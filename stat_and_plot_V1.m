    
    
    function statistic = stat_and_plot_V1( A, A_h, dM, beta, N_e, dim, M_slice, dt )
        
        N1 = size( A, 1 );    % N1 is the number of rows of the data storage matrix A
        N2 = size( A, 2 );    % N2 is the number of columns of the data storage matrix A
        M = N1 * N2;    % M is the total sample size of all the estimators for the Monte Carlo method
        M_power_max = round( log( M ) / log( 2 ) );
        M_power_min = 6;
        log2Mx = ( 6 : dM : M_power_max )';
        M_values = round( 2.^log2Mx );    % Mx_values are discrete points corresponding to the x-axis of the data points in the figure
        Num_M_points = length( M_values );    % N_Mx is the number of points to be plotted in the figure
        
        data_store_vec = zeros( M, 2 );
        for n1 = 1 : 1 : N1
            for n2 = 1 : 1 : N2
                m = ( n1 - 1 ) * N2 + n2;
                data_store_vec( m, 1 ) = A( n1, n2 );
                data_store_vec( m, 2 ) = A_h( n1, n2 );
                % data_store_vec( m, 2 ) = -A_h( n1, n2 ) / A( n1, n2 );
            end
        end
        
        Tr_path_integ_Z = zeros( Num_M_points, 1 );
        h_path_integ = zeros( Num_M_points, 1 );
        h_bar_std = zeros( Num_M_points, 1 );
        Z_bar_std = zeros( Num_M_points, 1 );
        
        Z_h_bar_cov = zeros( Num_M_points, 3 );
        h_bar_std_with_cov = zeros( Num_M_points, 1 );
        
        h_bar_moments = zeros( Num_M_points, 3 );
        Z_bar_moments = zeros( Num_M_points, 3 );
        
        for i_b = 1 : 1 : Num_M_points    % compute the mean of the MC integration samples
            M_ib = M_values( i_b, 1 );
            Tr_path_integ_Z( i_b, 1 ) = mean( data_store_vec( 1 : M_ib, 1 ), 1 );
            h_path_integ( i_b, 1 ) = -mean( data_store_vec( 1 : M_ib, 2 ), 1 ) / Tr_path_integ_Z( i_b, 1 );
            Z_bar_std( i_b, 1 ) = std( data_store_vec( 1 : M_ib, 1 ) );
            h_bar_std( i_b, 1 ) = std( ( data_store_vec( 1 : M_ib, 2 ) / Tr_path_integ_Z( i_b, 1 ) ) );
            
            Covariance_Z_h = cov( data_store_vec( 1 : M_ib, 1 ), data_store_vec( 1 : M_ib, 2 ) );
            Z_h_bar_cov( i_b, 1 ) = Covariance_Z_h( 1, 1 );
            Z_h_bar_cov( i_b, 2 ) = Covariance_Z_h( 2, 2 );
            Z_h_bar_cov( i_b, 3 ) = Covariance_Z_h( 1, 2 );
            
            Z_bar_moments( i_b, 1 ) = moment( data_store_vec( 1 : M_ib, 1 ), 2, 1 );
            Z_bar_moments( i_b, 2 ) = moment( data_store_vec( 1 : M_ib, 1 ), 3, 1 );
            Z_bar_moments( i_b, 3 ) = moment( data_store_vec( 1 : M_ib, 1 ), 4, 1 );
            h_bar_moments( i_b, 1 ) = moment( ( data_store_vec( 1 : M_ib, 2 ) / Tr_path_integ_Z( i_b, 1 ) ), 2, 1 );
            h_bar_moments( i_b, 2 ) = moment( ( data_store_vec( 1 : M_ib, 2 ) / Tr_path_integ_Z( i_b, 1 ) ), 3, 1 );
            h_bar_moments( i_b, 3 ) = moment( ( data_store_vec( 1 : M_ib, 2 ) / Tr_path_integ_Z( i_b, 1 ) ), 4, 1 );
        end
        
        % Two levels of for-loops using batches of samples
        Mxmax = log( N1 ) / log( 2 );
        N_estim = N2;
        dMx = 0.5;
        log2Mx = ( 6 : dMx : Mxmax )';
        Mx_values = round( 2.^log2Mx );    % Mx_values are discrete points corresponding to the x-axis of the data points in the figure

        N_Mx = length( Mx_values );    % N_Mx is the number of points to be plotted in the figure
        Mx_max = max( Mx_values );    % Mx_max is the maximum number of samples 
        % Allocate memory
        Tr_path_integ_Z_batch = zeros( N_Mx, N_estim );
        h_path_integ_batch = zeros( N_Mx, N_estim );
        for i_b = 1 : 1 : N_Mx    % compute the mean of the MC integration samples
            M_x = Mx_values( i_b, 1 );
            Tr_path_integ_Z_batch( i_b, : ) = mean( A( 1 : M_x, : ), 1 );
            h_path_integ_batch( i_b, : ) = -mean( A_h( 1 : M_x, : ), 1 ) ./ Tr_path_integ_Z_batch( i_b, : );
        end
        
        % Bootstrap on the samples of independent estimators of mean-field h, using mean function on the Bootstrap sample sets
        Num_bootstrap_sample = 400;
        [ h_bootstrap, bsamp ] = bootstrp( Num_bootstrap_sample, @(x) [ mean( x ) ] , h_path_integ_batch( N_Mx, : ) );
        h_bootstrap_CI = bootci( Num_bootstrap_sample, @(x) [ mean( x ) ] , h_path_integ_batch( N_Mx, : ) );
        
        
        % Plotting and outputting of data
        h_ref = h_eigval_Z_dq( N_e, beta, dim );
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        fprintf( 'Method 1, using all the samples of each independent Brownian bridges:\n' );
        fprintf( 'Inverse temperature beta = %.2f, %d electrons in dimension %d, M = %d time slices with Delta-t = %.3e, total sample size M_sp = %d \n', beta, N_e, dim, M_slice, dt, M );
        
        h_bar_Mmax = h_path_integ( Num_M_points, 1 );
        h_bar_std_Mmax = h_bar_std( Num_M_points, 1 );
        Z_bar_Mmax = Tr_path_integ_Z( Num_M_points, 1 );
        Z_bar_std_Mmax = Z_bar_std( Num_M_points, 1 );
        
        dZ_Z_cov_Mmax = Z_h_bar_cov( Num_M_points, 3 );
        
        h_bar_abs_CI = 1.96 * h_bar_std_Mmax / sqrt( M );
        h_bar_rel_CI = 1.96 * h_bar_std_Mmax / sqrt( M ) / h_bar_Mmax;
        % Considering the covariance between the numerator A and the denominator B, with A = dZ_dbeta, B = -Z.
        % The third term is -2 * Cov( A, B ) / ( A * B )
        h_bar_rel_std_with_cov = sqrt( ( h_bar_std_Mmax / h_bar_Mmax )^2 + ( Z_bar_std_Mmax / Z_bar_Mmax )^2 - 2 * dZ_Z_cov_Mmax / ( -Z_bar_Mmax * h_bar_Mmax * Z_bar_Mmax ) );
        h_bar_abs_CI_cov = 1.96 * h_bar_rel_std_with_cov * h_bar_Mmax / sqrt( M );
        h_bar_rel_CI_cov = 1.96 * h_bar_rel_std_with_cov / sqrt( M );
        
        h_rel_error = abs( h_bar_Mmax - h_ref ) / h_ref;
        h_re_CI_one_side = h_bar_std_Mmax / sqrt( M ) * 1.96 / h_ref;
        h_re_CI_one_side_cov = h_bar_rel_std_with_cov * h_bar_Mmax / sqrt( M ) * 1.96 / h_ref;
        
        fprintf( 'Reference value of mean-field energy h_ref = %6.5e\n', h_ref );
        fprintf( 'Using the determinant formula, applying method 1 without considering covariance, approximate value of the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax, h_bar_abs_CI, h_bar_rel_CI );
        fprintf( 'h_bar has relative error %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error, h_re_CI_one_side );
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        fprintf( 'Using the determinant formula, applying method 1 with covariance between numerator and denominator, approximate value of the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax, h_bar_abs_CI_cov, h_bar_rel_CI_cov );
        fprintf( 'h_bar has relative error %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error, h_re_CI_one_side_cov );
        
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        
        % Output data: beta, N_e, dim, M_time_slice, dt, h_ref, Mxmax, N_estim, h_bar_Mmax, h_bar_abs_CI, h_rel_diff_with_ref, h_re_CI_one_side
		if ( ~exist( 'Table_info', 'dir' ) )
			mkdir( 'Table_info' );
		end
        if ( exist( 'Table_info/Pure_HO_Approx_Z_and_h.txt', 'file' ) == 0 )
            fileID = fopen( "Table_info/Pure_HO_Approx_Z_and_h.txt", 'a' );
            fprintf( fileID, 'beta,  N_e,  dim,  M_time_slice,  dt,  h_ref,  Mxmax,  N_estim,  h_bar_Mmax,  h_bar_abs_CI,  h_rel_diff_with_ref,  h_re_CI_one_side\n' );
            fclose( fileID );
        end        
        fileID = fopen( "Table_info/Pure_HO_Approx_Z_and_h.txt", 'a' );
        fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
        fprintf( fileID, 'Method 1, without covariance, mean-field h:           %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, h_ref, Mxmax, N_estim, h_bar_Mmax, h_bar_abs_CI, h_rel_error, h_re_CI_one_side );
        fprintf( fileID, 'Method 1, with covariance, mean-field h:              %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, h_ref, Mxmax, N_estim, h_bar_Mmax, h_bar_abs_CI_cov, h_rel_error, h_re_CI_one_side_cov );
        fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
        fclose( fileID );
        
        
        fig_1 = figure;
        loglog( M_values( 1 : 1 : Num_M_points ), abs( h_path_integ( 1 : 1 : Num_M_points ) - h_ref ) / h_ref, 'o-' );
        hold on
        loglog( M_values, M_values.^(-0.5) );
        lgd = legend( 'Relative difference in approx. of h with varing M_{sp} value', 'reference line of M_{sp}^{-1/2}' );
        lgd.FontSize = 20;
        title( "n = " + N_e + " fermions in dim " + dim + ", computing h by taking derivative on the approx. with determinant formula, Coulomb + HO potential, \beta = " + beta + ", \Delta t = " + dt );
        xlabel( 'Number of sample points M_{sp} for MC integral' );
        ylabel( 'Rel-error of Monte Carlo integral' );
        ax = gca;
        ax.FontSize = 20;
        hold off
        
        fig_3 = figure;
        plot( M_values( 1 : 1 : Num_M_points ), h_path_integ( 1 : 1 : Num_M_points ), 'o-' );
        hold on
        plot( M_values( 1 : 1 : Num_M_points ), h_bar_moments( 1 : 1 : Num_M_points, 1 ), 'o-' );
        plot( M_values( 1 : 1 : Num_M_points ), h_bar_moments( 1 : 1 : Num_M_points, 2 ), 'o-' );
        plot( M_values( 1 : 1 : Num_M_points ), h_bar_moments( 1 : 1 : Num_M_points, 3 ), 'o-' );
        lgd = legend( 'mean value', '2nd moment', '3rd moment', '4th moment' );
        lgd.FontSize = 20;
        title( "n = " + N_e + " fermions in dim " + dim + ", approx. of h, Coulomb + HO potential, \beta = " + beta + ", \Delta t = " + dt );
        xlabel( 'Number of sample points M_{sp} for MC integral' );
        ylabel( 'Mean and moments of Monte Carlo integral estimators for h' );
        ax = gca;
        ax.FontSize = 20;
        hold off
        
        
        Z_ref = Z_partition_a( N_e, beta, dim );
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        Z_bar_Mmax = Tr_path_integ_Z( Num_M_points, 1 );
        Z_bar_std_Mmax = Z_bar_std( Num_M_points, 1 );
        Z_bar_abs_CI = 1.96 * Z_bar_std_Mmax / sqrt( M );
        Z_bar_rel_CI = 1.96 * Z_bar_std_Mmax / sqrt( M ) / Z_bar_Mmax;
        Z_rel_error = abs( Z_bar_Mmax - Z_ref ) / Z_ref;
        Z_re_CI_one_side = Z_bar_std_Mmax / sqrt( M ) * 1.96 / Z_ref;
        fprintf( 'Reference value of partition function Z_ref = %6.5e\n', Z_ref );
        fprintf( 'Using the determinant formula to approximate the partition function Z_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', Z_bar_Mmax, Z_bar_abs_CI, Z_bar_rel_CI );
        fprintf( 'Z_bar has relative error %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', Z_rel_error, Z_re_CI_one_side );
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        fileID = fopen( "Table_info/Pure_HO_Approx_Z_and_h.txt", 'a');
        fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
        fprintf( fileID, 'Method 1, partition function Z:                       %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, Z_ref, Mxmax, N_estim, Z_bar_Mmax, Z_bar_abs_CI, Z_rel_error, Z_re_CI_one_side );
        fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
        fclose( fileID );

        fig_2 = figure;
        loglog( M_values( 1 : 1 : Num_M_points ), abs( Tr_path_integ_Z( 1 : 1 : Num_M_points ) - Z_ref ) / Z_ref, 'o-' );
        hold on
        loglog( M_values, M_values.^(-0.5) );
        lgd = legend( 'Relative difference in approx. of Z with varing M_{sp} value', 'reference line of M_{sp}^{-1/2}' );
        lgd.FontSize = 20;
        title( "n = " + N_e + " fermions in dim " + dim + ", computing Z with determinant formula, pure HO potential, \beta = " + beta + ", \Delta t = " + dt );
        xlabel( 'Number of sample points M_{sp} for MC integral' );
        ylabel( 'Rel-error of Monte Carlo integral' );
        ax = gca;
        ax.FontSize = 20;
        hold off
        
        fig_4 = figure;
        plot( M_values( 1 : 1 : Num_M_points ), Tr_path_integ_Z( 1 : 1 : Num_M_points ), 'o-' );
        hold on
        plot( M_values( 1 : 1 : Num_M_points ), Z_bar_moments( 1 : 1 : Num_M_points, 1 ), 'o-' );
        plot( M_values( 1 : 1 : Num_M_points ), Z_bar_moments( 1 : 1 : Num_M_points, 2 ), 'o-' );
        plot( M_values( 1 : 1 : Num_M_points ), Z_bar_moments( 1 : 1 : Num_M_points, 3 ), 'o-' );
        lgd = legend( 'mean value', '2nd moment', '3rd moment', '4th moment' );
        lgd.FontSize = 20;
        title( "n = " + N_e + " fermions in dim " + dim + ", approx. of Z, Coulomb + HO potential, \beta = " + beta + ", \Delta t = " + dt );
        xlabel( 'Number of sample points M_{sp} for MC integral' );
        ylabel( 'Mean and moments of Monte Carlo integral estimators for Z' );
        ax = gca;
        ax.FontSize = 20;
        hold off
        
		if ( ~exist( 'Figures_store', 'dir' ) )
			mkdir( 'Figures_store' );
		end
        saveas( fig_1, [ pwd sprintf( '/Figures_store/Pure_HO_Approx_h_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
        close( fig_1 );
        saveas( fig_3, [ pwd sprintf( '/Figures_store/Pure_HO_Moments_h_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
        close( fig_3 );
        
        saveas( fig_2, [ pwd sprintf( '/Figures_store/Pure_HO_Approx_Z_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
        close( fig_2 );
        saveas( fig_4, [ pwd sprintf( '/Figures_store/Pure_HO_Moments_Z_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
        close( fig_4 );
        
        
        h_bar_approx = mean( h_path_integ_batch, 2 );
        h_bar_std_batch = std( h_path_integ_batch, 0, 2 );
        h_bar_Mmax_batch = h_bar_approx( N_Mx, 1 );
        h_bar_std_Mmax_batch = h_bar_std_batch( N_Mx, 1 );
        h_bar_abs_CI_batch = 1.96 * h_bar_std_Mmax_batch / sqrt( N_estim );
        h_bar_rel_CI_batch = 1.96 * h_bar_std_Mmax_batch / sqrt( N_estim ) / h_bar_Mmax_batch;
        h_rel_error_batch = abs( h_bar_Mmax_batch - h_ref ) / h_ref;
        h_re_CI_one_side_batch = h_bar_std_Mmax_batch / sqrt( N_estim ) * 1.96 / h_ref;
        
        h_mean_bootstrap = mean( h_bootstrap );
        h_bar_bootstrap_abs_CI_built_in = abs( h_bootstrap_CI( 2, 1 ) - h_bootstrap_CI( 1, 1 ) ) / 2;
        h_bar_bootstrap_rel_CI_built_in = h_bar_bootstrap_abs_CI_built_in / h_mean_bootstrap;
        h_rel_error_bootstrap = abs( h_mean_bootstrap - h_ref ) / h_ref;
        
        
        fprintf( 'Reference value of mean-field energy h_ref = %6.5e\n', h_ref );
        fprintf( 'Using the determinant formula, applying method 2 without bootstrapping, approximate value of the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax_batch, h_bar_abs_CI_batch, h_bar_rel_CI_batch );
        fprintf( 'h_bar has relative error %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error_batch, h_re_CI_one_side_batch );
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        fprintf( 'Using the determinant formula, applying the method 2 with bootstraping to approximate the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e \n', h_mean_bootstrap, h_bar_bootstrap_abs_CI_built_in );
        fprintf( 'h_mean_bootstrap has relative error %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error_bootstrap, h_bar_bootstrap_rel_CI_built_in );
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        
        fprintf( '—————————————————————————————————————————————————————————————————————\n' );
        fileID = fopen( "Table_info/Pure_HO_Approx_Z_and_h.txt", 'a');
        fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
        fprintf( fileID, 'Method 2, without bootstraping, mean-field h:         %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, h_ref, Mxmax, N_estim, h_bar_Mmax_batch, h_bar_abs_CI_batch, h_rel_error_batch, h_re_CI_one_side_batch );
        fprintf( fileID, 'Method 2, with bootstraping, mean-field h:            %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, h_ref, Mxmax, N_estim, h_mean_bootstrap, h_bar_bootstrap_abs_CI_built_in, h_rel_error_bootstrap, h_bar_bootstrap_rel_CI_built_in );
        fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
        fclose( fileID );
        
        
        
        
        
        
        
        
    
        statistic = 1;
    end