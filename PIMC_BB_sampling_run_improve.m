
function exitcode = PIMC_BB_sampling_run_improve( filename0, lulog0, datadir0, ...
    dim0, Ne0, beta0, lambda0, dt0, Mxmax0, dMx0, Nrep0, Ncore0, TestCase0 )
%
% Stand-alone function for path integral MD sampling.
%
%
exitcode = 0;
%
savechkp = 1;     % save checkpoints
%
% Set Parallel pool
Ncore = str2num(Ncore0); 
myPool = parpool( Ncore );
%
% Reading cmd line parameters
if nargin == 0
    fprintf( 1, 'Usage: \n' );
    fprintf( 1, '    exitcode  = sampling_MP_ULD(outfile,Ntr, Tfinal,at,bt,nt)\n' );
    fprintf( 1, '    filename  = name of output files \n' );
    fprintf( 1, '    lulog     = output messages: 0... stdout, 1... fileprefix.log \n' );
    fprintf( 1, '    datadir   = directory with data files \n' );
    fprintf( 1, '    dim0      = dimensionality of the space. \n');
    fprintf( 1, '    Ne0       = Number of electrons. \n');
    fprintf( 1, '    beta0     = beta=inverse temp. \n');
    fprintf( 1, '    lambda0   = lambda parameter for characterising the Coulomb strength. \n');
    fprintf( 1, '    dt        = time step in Trotter expansion \n' );
    fprintf( 1, '    Mxmax0    = log2 for max Mx \n' );
    fprintf( 1, '    dMx0      = log2 delta Mx \n' );
    fprintf( 1, '    Nrep0   = number of independent runs \n' );

    fprintf( 1, '\n' );
    exitcode = -1;
    return;
end
if isdeployed == 1
    fprintf(1, '===> [PIMC_BB_samplig_run] is deployed\n');
end

fprintf(1, '==> PARALLEL POOL: NumWorkers: %d\n',myPool.NumWorkers);

numprocs = 2;

dirname = datadir0;
outfile = [filename0,'.mat'];
lulog   = str2num(lulog0);
dim     = str2num(dim0);
N_e     = str2num(Ne0);
beta    = str2num(beta0);        % beta in the manuscript
lambda  = str2num(lambda0);
dt      = str2num(dt0);          % dt is the \Delta t value of the manuscript
Mxmax   = str2num(Mxmax0);
dMx     = str2num(dMx0);
N_estim = 2^str2num(Nrep0);       % N_estim is the number of realisations 
TestCase = TestCase0;
cpuflag = 1;
%
% system parameters hard-wired 
d = dim * N_e;    
M_sp = 1;
M_slice = round( beta / dt );    % M_slice is the number of beads

log2Mx = ( 6 : dMx : Mxmax )';
Mx_values = round( 2.^log2Mx );    % Mx_values are discrete points corresponding to the x-axis of the data points in the figure

N_Mx = length( Mx_values );    % N_Mx is the number of points to be plotted in the figure
Mx_max = max( Mx_values );    % Mx_max is the maximum number of samples generated with overdamped Langevin method

% Allocate memory
MC_integrand_Z_sp_store = zeros( Mx_max, N_estim );
Tr_path_integ_Z = zeros( N_Mx, N_estim );  
MC_integrand_h_sp_store = zeros( Mx_max, N_estim );
h_path_integ = zeros( N_Mx, N_estim );

tini_sample = tic;
ticBytes(gcp);

if( TestCase == "V1" )
    lambda = 0;
	% Loop over independent realization
	parfor ( ell = 1 : N_estim )
		
		task = getCurrentTask();
		if( mod( ell, 100 ) == 0 )
			fprintf(1, '--> [BGN: %d] replica of sample set %d\n',task.ID,ell);
		end
		
		two_poss_dens_rv = rand( Mx_max, 1 );
		% two_poss_dens_rv is a uniform r.v. in [0,1], is it is smaller than 1/2, we use Normal distribution with variance \beta for the sampling of initial points x0;
		% Otherwise if two_poss_dens_rv is larger than 1/2, we generate initial points x0 under the Normal distribution with variance 1/beta.    
		
		for p = 1 : 1 : Mx_max
			
			[ x_0, p_x_density ] = sample_x0( beta, dim, N_e );
			
			sum_over_sp = 0;    % sum_over_sp stores the sum of the weighted path integral values exp(-\int V(x_t) dt) * weight(x) for each sample path
			sum_over_sp_dZ = 0;
			
			for j = 1 : 1 : M_sp
				W_incre_store = randn( d, M_slice );
				% preparing Brownian increments of each Brownian bridge
				% rng( 127 )
				W_incre = sqrt( dt ) * W_incre_store;
				W_t = zeros( d, M_slice + 1 );
				% Generating the Brownian paths
				sum_temp_Brownian_path = zeros( d, 1 );
				for i = 1 : 1 : M_slice
					sum_temp_Brownian_path = sum_temp_Brownian_path + W_incre( :, i );
					W_t( :, i + 1 ) = sum_temp_Brownian_path;
				end
				B_t = zeros( d, M_slice + 1 );
				W_T = W_t( :, M_slice + 1 );
				for i = 1 : 1 : M_slice
					W_t_i = W_t( :, i );
					B_t_i = W_t_i - ( ( i - 1 ) / M_slice ) * W_T;
					B_t( :, i ) = B_t_i;
				end
				B_t( :, M_slice + 1 ) = B_t( :, 1 );
				
				[ W_mat, dW_dbeta_mat ] = W_matrices_mod_V1( x_0, B_t, beta, dt, N_e, M_slice, dim, lambda );
				det_W = det( W_mat );
				adjugate_W = adjoint( W_mat );
				trace_dW = sum( diag( adjugate_W * dW_dbeta_mat ) );
				summand_over_sp_dZ = trace_dW - det_W * dim * N_e / ( 2 * beta );
				sum_over_sp = sum_over_sp + det_W;    % compute the sum of weighted path integrals, at inverse temperature beta
				sum_over_sp_dZ = sum_over_sp_dZ + summand_over_sp_dZ;    % compute the sum of weighted path integrals, at \beta + \Delta\beta
				
			end
			mean_over_sp = sum_over_sp / M_sp * ( 2 * pi * beta )^( -( dim * N_e ) / 2 ) / factorial( N_e );
			g_x_p = mean_over_sp / p_x_density;    % normalize the samples with the corresponding probability density of the initial point x0
			MC_integrand_Z_sp_store( p, ell ) = g_x_p;    % store the samples for the Monte-Carlo integration
			mean_over_sp_dZ = sum_over_sp_dZ / M_sp * ( 2 * pi * beta )^( -( dim * N_e ) / 2 ) / factorial( N_e );
			g_x_p_dZ = mean_over_sp_dZ / p_x_density;
			MC_integrand_h_sp_store( p, ell ) = g_x_p_dZ;    % store the samples for the Monte-Carlo integration
					
		end
		
		if( mod( ell, 100 ) == 0 )
			fprintf(1,'--> [END: %d] replica of sample set \n',task.ID);
		end
	end
elseif( TestCase == "V2" )
	% Loop over independent realization
	parfor ( ell = 1 : N_estim )
		
		task = getCurrentTask();
		if( mod( ell, 100 ) == 0 )
			fprintf(1, '--> [BGN: %d] replica of sample set %d\n',task.ID,ell);
		end
		
		two_poss_dens_rv = rand( Mx_max, 1 );
		% two_poss_dens_rv is a uniform r.v. in [0,1], is it is smaller than 1/2, we use Normal distribution with variance \beta for the sampling of initial points x0;
		% Otherwise if two_poss_dens_rv is larger than 1/2, we generate initial points x0 under the Normal distribution with variance 1/beta.    
		
		for p = 1 : 1 : Mx_max
			
			[ x_0, p_x_density ] = sample_x0( beta, dim, N_e );
			
			sum_over_sp = 0;    % sum_over_sp stores the sum of the weighted path integral values exp(-\int V(x_t) dt) * weight(x) for each sample path
			sum_over_sp_dZ = 0;
			
			for j = 1 : 1 : M_sp
				W_incre_store = randn( d, M_slice );
				% preparing Brownian increments of each Brownian bridge
				% rng( 127 )
				W_incre = sqrt( dt ) * W_incre_store;
				W_t = zeros( d, M_slice + 1 );
				% Generating the Brownian paths
				sum_temp_Brownian_path = zeros( d, 1 );
				for i = 1 : 1 : M_slice
					sum_temp_Brownian_path = sum_temp_Brownian_path + W_incre( :, i );
					W_t( :, i + 1 ) = sum_temp_Brownian_path;
				end
				B_t = zeros( d, M_slice + 1 );
				W_T = W_t( :, M_slice + 1 );
				for i = 1 : 1 : M_slice
					W_t_i = W_t( :, i );
					B_t_i = W_t_i - ( ( i - 1 ) / M_slice ) * W_T;
					B_t( :, i ) = B_t_i;
				end
				B_t( :, M_slice + 1 ) = B_t( :, 1 );
				
				[ W_mat, dW_dbeta_mat ] = W_matrices_mod_V2( x_0, B_t, beta, dt, N_e, M_slice, dim, lambda );
				det_W = det( W_mat );
				adjugate_W = adjoint( W_mat );
				trace_dW = sum( diag( adjugate_W * dW_dbeta_mat ) );
				summand_over_sp_dZ = trace_dW - det_W * dim * N_e / ( 2 * beta );
				sum_over_sp = sum_over_sp + det_W;    % compute the sum of weighted path integrals, at inverse temperature beta
				sum_over_sp_dZ = sum_over_sp_dZ + summand_over_sp_dZ;    % compute the sum of weighted path integrals, at \beta + \Delta\beta
				
			end
			mean_over_sp = sum_over_sp / M_sp * ( 2 * pi * beta )^( -( dim * N_e ) / 2 ) / factorial( N_e );
			g_x_p = mean_over_sp / p_x_density;    % normalize the samples with the corresponding probability density of the initial point x0
			MC_integrand_Z_sp_store( p, ell ) = g_x_p;    % store the samples for the Monte-Carlo integration
			mean_over_sp_dZ = sum_over_sp_dZ / M_sp * ( 2 * pi * beta )^( -( dim * N_e ) / 2 ) / factorial( N_e );
			g_x_p_dZ = mean_over_sp_dZ / p_x_density;
			MC_integrand_h_sp_store( p, ell ) = g_x_p_dZ;    % store the samples for the Monte-Carlo integration
					
		end
		
		if( mod( ell, 100 ) == 0 )
			fprintf(1,'--> [END: %d] replica of sample set \n',task.ID);
		end
	end
end

% Implement the statistics and generate the figures
if( TestCase == "V1" )
	stat_plot = stat_and_plot_V1( MC_integrand_Z_sp_store, MC_integrand_h_sp_store, dMx, beta, N_e, dim, M_slice, dt )
elseif( TestCase == "V2" )
    stat_plot = stat_and_plot_V2( MC_integrand_Z_sp_store, MC_integrand_h_sp_store, dMx, beta, N_e, dim, M_slice, dt )
end

sample_cpu = toc( tini_sample );
if cpuflag
    fprintf(1,'CPU per %d samples: %e [s]\n',N_estim,sample_cpu);
    [hf,mf,sf] = hms(seconds(sample_cpu));
    fprintf(1,'Estimated CPU: %d:%d:%4.2f\n',hf,mf,sf);
    cpuflag=0;
end

tocBytes(gcp)

delete(myPool);
fprintf(1,'[END] PIMC_BB_sample_run completed: data saved to %s\n',...
    outfile);
    

end

    
    
    
    
    
    
    