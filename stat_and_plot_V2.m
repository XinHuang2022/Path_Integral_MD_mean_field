    
    
    function statistic = stat_and_plot_V2( A, A_h, dM, beta, N_e, dim, M_slice, dt )
        
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
        
        
        % Reference value for dim = 2, beta = 1, N_e varing from 3 to 10
        if( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 3 ) )
            h_ref = 8.719;    % Dornheim's reference result for N_e = 3, dim = 2, beta = 1, h = 8.719(3)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 4 ) )
            h_ref = 12.903;    % Dornheim's reference result for N_e = 4, dim = 2, beta = 1, h = 12.903(7)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 5 ) )
            h_ref = 17.66;    % Dornheim's reference result for N_e = 5, dim = 2, beta = 1, h = 17.66(2)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 22.82;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 1, h = 22.82(5)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 7 ) )
            h_ref = 28.7;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 1, h = 28.7(2)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 8 ) )
            h_ref = 34.5;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 1, h = 34.5(5)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 9 ) )
            h_ref = 40;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 1, h = 40(2)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 2 ) && ( N_e == 10 ) )
            h_ref = 49;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 1, h = 49(3)
            Reference_exist = 1;
        % Reference value for dim = 2, beta = 0.3, N_e varying from 4 to 20
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 4 ) )
            h_ref = 29.374;    % Dornheim's reference result for N_e = 4, dim = 2, beta = 0.3, h = 29.374(8)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 46.45;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 0.3, h = 46.45(1)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 8 ) )
            h_ref = 64.97;    % Dornheim's reference result for N_e = 8, dim = 2, beta = 0.3, h = 64.97(2)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 10 ) )
            h_ref = 84.92;    % Dornheim's reference result for N_e = 10, dim = 2, beta = 0.3, h = 84.92(4)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 12 ) )
            h_ref = 106.20;    % Dornheim's reference result for N_e = 12, dim = 2, beta = 0.3, h = 106.20(6)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 16 ) )
            h_ref = 152.4;    % Dornheim's reference result for N_e = 16, dim = 2, beta = 0.3, h = 152.4(3)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 18 ) )
            h_ref = 179.1;    % Dornheim's reference result for N_e = 18, dim = 2, beta = 0.3, h = 179.1(6)
            Reference_exist = 1;
        elseif( ( beta == 0.3 ) && ( dim == 2 ) && ( N_e == 20 ) )
            h_ref = 203.0;    % Dornheim's reference result for N_e = 20, dim = 2, beta = 0.3, h = 203(1)
            Reference_exist = 1;
        % Reference value for dim = 2, N_e = 6, beta varying from 0.3 to 2
        elseif( ( beta == 0.5 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 32.16;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 0.5, h = 32.16(1)
            Reference_exist = 1;
        elseif( ( beta == 0.6 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 28.86;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 0.6, h = 28.86(1)
            Reference_exist = 1;
        elseif( ( beta == 0.8 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 24.98;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 0.8, h = 24.98(3)
            Reference_exist = 1;
        elseif( ( beta == 1.5 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 20.4;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 1.5, h = 20.4(5)
            Reference_exist = 1;
        elseif( ( beta == 2 ) && ( dim == 2 ) && ( N_e == 6 ) )
            h_ref = 17;    % Dornheim's reference result for N_e = 6, dim = 2, beta = 2, h = 17(1)
            Reference_exist = 1;
        % Reference value for dim = 3, N_e = 6, beta varying from 0.3 to 2
        elseif( ( beta == 0.3 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 64.09;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 0.3, h = 64.09(2)
            Reference_exist = 1;
        elseif( ( beta == 0.5 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 41.66;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 0.5, h = 41.66(1)
            Reference_exist = 1;
        elseif( ( beta == 0.6 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 36.341;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 0.6, h = 36.341(9)
            Reference_exist = 1;
        elseif( ( beta == 0.8 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 30.114;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 0.8, h = 30.114(9)
            Reference_exist = 1;
        elseif( ( beta == 1 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 26.692;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 1, h = 26.692(9)
            Reference_exist = 1;
        elseif( ( beta == 1.5 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 22.63;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 1.5, h = 22.63(7)
            Reference_exist = 1;
        elseif( ( beta == 2 ) && ( dim == 3 ) && ( N_e == 6 ) )
            h_ref = 22.1;    % Dornheim's reference result for N_e = 6, dim = 3, beta = 1, h = 22.1(5)
            Reference_exist = 1;
        % Reference value computed from taking the derivative w.r.t. beta on the tensor expression
        elseif( ( beta == 1 ) && ( dim == 3 ) && ( N_e == 3 ) )
            h_ref = 11.355;    % Reference value obtained from the tensor expression
            Reference_exist = 1;
        elseif( ( beta == 1.5 ) && ( dim == 3 ) && ( N_e == 3 ) )
            h_ref = 9.157;    % Reference value obtained from the tensor expression
            Reference_exist = 1;
        elseif( ( beta == 2 ) && ( dim == 3 ) && ( N_e == 3 ) )
            h_ref = 8.282;    % Reference value obtained from the tensor expression
            Reference_exist = 1;
        else
            Reference_exist = 0;
            
        end
        

        if( Reference_exist == 0 )
            fprintf( '—————————————————————————————————————————————————————————————————————————————\n' );
            fprintf( 'Inverse temperature beta = %.2f, %d electrons in dimension %d, M = %d time slices and Delta-t = %.3e, M_sp = %d \n', beta, N_e, dim, M_slice, dt, M );
            h_bar_Mmax = h_path_integ( Num_M_points, 1 );
            h_bar_std_Mmax = h_bar_std( Num_M_points, 1 );
            h_bar_abs_CI = 1.96 * h_bar_std_Mmax / sqrt( M );
            h_bar_rel_CI = 1.96 * h_bar_std_Mmax / sqrt( M ) / h_bar_Mmax;
            fprintf( 'Using the determinant formula to approximate the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax, h_bar_abs_CI, h_bar_rel_CI );
            fprintf( '—————————————————————————————————————————————————————————————————————————————\n' );
            fileID = fopen( "Table_info/h_store_beta=" + beta + "_dim=" + dim + "_n=" + N_e + "_dt=" + dt + "_M_sp=" + M + ".txt", 'w');
            fprintf( fileID, '—————————————————————————————————————————————————————————————————————————————\n' );
            fprintf( fileID, 'Inverse temperature beta = %.2f, %d electrons in dimension %d, M = %d time slices and Delta-t = %.3e \n', beta, N_e, dim, M_slice, dt );
            fprintf( fileID, 'Using the determinant formula to approximate the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax, h_bar_abs_CI, h_bar_rel_CI );
            fprintf( fileID, '—————————————————————————————————————————————————————————————————————————————\n' );
            fclose( fileID );
            
            fig_1 = figure;
            loglog( M_values( 1 : 1 : Num_M_points - 1 ), abs( h_path_integ( 1 : 1 : Num_M_points - 1 ) - h_bar_Mmax ) / h_bar_Mmax, 'o-' );
            hold on
            loglog( M_values, M_values.^(-0.5) );
            lgd = legend( 'Relative difference in approx. of h with varing M_x value', 'reference line of M_x^{-1/2}' );
            lgd.FontSize = 20;
            title( "n = " + N_e + " fermions in dim " + dim + ", computing h by taking derivative on the approx. with determinant formula, Coulomb + HO potential, \beta = " + beta + ", \Delta t = " + dt );
            xlabel( 'Number of sample points M_x for MC integral' );
            ylabel( 'Rel-error of Monte Carlo integral' );
            ax = gca;
            ax.FontSize = 20;
            hold off
            
			if ( ~exist( 'Figures_store', 'dir' ) )
				mkdir( 'Figures_store' );
			end
            saveas( fig_1, [ pwd sprintf( '/Figures_store/Approx_h_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_no_ref.fig', beta, N_e, dim, dt, M ) ] );
            close( fig_1 );
        end


        if( Reference_exist == 1 )
        
            % Plotting and outputting of data
            fprintf( '—————————————————————————————————————————————————————————————————————\n' );
            fprintf( 'Method 1, using all the samples of each independent Brownian bridges:\n' );
            fprintf( 'Inverse temperature beta = %.2f, %d electrons in dimension %d, M = %d time slices with Delta-t = %.3e, total sample size M_sp = %d \n', beta, N_e, dim, M_slice, dt, M );
            % fprintf( 'Each estimator employing 2 to the power %d M_x0 sample points, using %d independent estimators in total, to approximately compute the statistical error in h \n', Mxmax, N_estim );
            
            h_bar_Mmax = h_path_integ( Num_M_points, 1 );
            h_bar_std_Mmax = h_bar_std( Num_M_points, 1 );
            Z_bar_Mmax = Tr_path_integ_Z( Num_M_points, 1 );
            Z_bar_std_Mmax = Z_bar_std( Num_M_points, 1 );
            
            dZ_Z_cov_Mmax = Z_h_bar_cov( Num_M_points, 3 );
            
            h_bar_abs_CI = 1.96 * h_bar_std_Mmax / sqrt( M );
            h_bar_rel_CI = 1.96 * h_bar_std_Mmax / sqrt( M ) / h_bar_Mmax;
            h_bar_rel_std_with_cov = sqrt( ( h_bar_std_Mmax / h_bar_Mmax )^2 + ( Z_bar_std_Mmax / Z_bar_Mmax )^2 - 2 * dZ_Z_cov_Mmax / ( -Z_bar_Mmax * h_bar_Mmax * Z_bar_Mmax ) );
            h_bar_abs_CI_cov = 1.96 * h_bar_rel_std_with_cov * h_bar_Mmax / sqrt( M );
            h_bar_rel_CI_cov = 1.96 * h_bar_rel_std_with_cov / sqrt( M );
            h_rel_error = abs( h_bar_Mmax - h_ref ) / h_ref;
            h_re_CI_one_side = h_bar_std_Mmax / sqrt( M ) * 1.96 / h_ref;
            h_re_CI_one_side_cov = h_bar_rel_std_with_cov * h_bar_Mmax / sqrt( M ) * 1.96 / h_ref;
            
            fprintf( 'Reference value of mean-field energy h_ref = %6.5e\n', h_ref );
            fprintf( 'Using the determinant formula, applying method 1 without considering covariance, approximate value of the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax, h_bar_abs_CI, h_bar_rel_CI );
            fprintf( 'h_bar has relative difference %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error, h_re_CI_one_side );
            fprintf( '—————————————————————————————————————————————————————————————————————\n' );
            fprintf( 'Using the determinant formula, applying method 1 with covariance between numerator and denominator, approximate value of the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e and relative CI %7.6e\n', h_bar_Mmax, h_bar_abs_CI_cov, h_bar_rel_CI_cov );
            fprintf( 'h_bar has relative difference %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error, h_re_CI_one_side_cov );
                
            fprintf( '—————————————————————————————————————————————————————————————————————\n' );
            
            % Output data: beta, N_e, dim, M_time_slice, dt, h_ref, Mxmax, N_estim, h_bar_Mmax, h_bar_abs_CI, h_rel_diff_with_ref, h_re_CI_one_side
            if ( ~exist( 'Table_info', 'dir' ) )
				mkdir( 'Table_info' );
			end
			if ( exist( 'Table_info/HO_and_Coulomb_Approx_Z_and_h.txt', 'file' ) == 0 )
                fileID = fopen( "Table_info/HO_and_Coulomb_Approx_Z_and_h.txt", 'a' );
                fprintf( fileID, 'beta,  N_e,  dim,  M_time_slice,  dt,  h_ref,  Mxmax,  N_estim,  h_bar_Mmax,  h_bar_abs_CI,  h_rel_diff_with_ref,  h_re_CI_one_side\n' );
                fclose( fileID );
            end
            fileID = fopen( "Table_info/HO_and_Coulomb_Approx_Z_and_h.txt", 'a');
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
			
            saveas( fig_1, [ pwd sprintf( '/Figures_store/HO_and_Coulomb_Approx_h_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
            close( fig_1 );
            saveas( fig_3, [ pwd sprintf( '/Figures_store/HO_and_Coulomb_Moments_h_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
            close( fig_3 );
            
            saveas( fig_4, [ pwd sprintf( '/Figures_store/HO_and_Coulomb_Moments_Z_beta=%.2f_n=%d_dim=%d_dt=%.3e_M_sp=%d_HO.fig', beta, N_e, dim, dt, M ) ] );
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
            fprintf( 'h_bar has relative difference %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error_batch, h_re_CI_one_side_batch );
            fprintf( '—————————————————————————————————————————————————————————————————————\n' );
            fprintf( 'Using the determinant formula, applying the method 2 with bootstraping to approximate the mean-field h_bar = %7.6e, with 0.95 absolute CI %7.6e \n', h_mean_bootstrap, h_bar_bootstrap_abs_CI_built_in );
            fprintf( 'h_mean_bootstrap has relative difference %7.6e, with 0.95 confidence interval plus / minus % 7.6e\n', h_rel_error_bootstrap, h_bar_bootstrap_rel_CI_built_in );
            fprintf( '—————————————————————————————————————————————————————————————————————\n' );
            
            fprintf( '—————————————————————————————————————————————————————————————————————\n' );
            fileID = fopen( "Table_info/HO_and_Coulomb_Approx_Z_and_h.txt", 'a');
            fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
            fprintf( fileID, 'Method 2, without bootstraping, mean-field h:         %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, h_ref, Mxmax, N_estim, h_bar_Mmax_batch, h_bar_abs_CI_batch, h_rel_error_batch, h_re_CI_one_side_batch );
            fprintf( fileID, 'Method 2, with bootstraping, mean-field h:            %.2f,   %d,   %d,   %d,   %.3e,   %6.5e,   2^%d,    %d,    %7.6e,    %7.6e,    %5.4e,    %5.4e \n', beta, N_e, dim, M_slice, dt, h_ref, Mxmax, N_estim, h_mean_bootstrap, h_bar_bootstrap_abs_CI_built_in, h_rel_error_bootstrap, h_bar_bootstrap_rel_CI_built_in );
            fprintf( fileID, '——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————\n' );
            fclose( fileID );
        
        end
        
        statistic = 1;
    end