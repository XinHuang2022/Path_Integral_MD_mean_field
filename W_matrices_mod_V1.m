    
    function [ W_mat, W_deriv_mat ] = W_matrices_mod_V1( x_0, B_t, beta, dt, N_e, M, dim, lambda )
    
        W_mat = zeros( N_e, N_e );
		W_deriv_mat = zeros( N_e, N_e );
        if( dim == 3 )
			x0_x = x_0( 1 : N_e, 1 );  x0_y = x_0( N_e + 1 : 2 * N_e, 1 );  x0_z = x_0( 2 * N_e + 1 : 3 * N_e, 1 );
			x0_vectors = [ x0_x, x0_y, x0_z ];
        elseif( dim == 2 )
            x0_x = x_0( 1 : N_e, 1 );  x0_y = x_0( N_e + 1 : 2 * N_e, 1 );
			x0_vectors = [ x0_x, x0_y ];
        elseif( dim == 1 )
            x0_x = x_0( 1 : N_e, 1 ); 
			x0_vectors = x0_x;
        end
        
        for k = 1 : 1 : N_e
            x0_k = x0_vectors( k, : );
            for ell = 1 : 1 : N_e
                x0_ell = x0_vectors( ell, : );
                x0_diff = x0_k - x0_ell;
                x0_diff_square = x0_diff * x0_diff';
                V_HO_sum = dt * ( V_func_HO( x0_k ) + V_func_HO( x0_ell ) ) / 2;
                for m = 2 : 1 : M
                    B_t_m = B_t( :, m );
                    B_k_t_m = BB_component( B_t_m, k, N_e, dim );
                    x_k_ell_t_m = B_k_t_m + ( 1 - ( m - 1 ) / M ) * x0_k + ( ( m - 1 ) / M ) * x0_ell;
                    V_HO_term = V_func_HO( x_k_ell_t_m );
                    V_HO_sum = V_HO_sum + dt * V_HO_term;
                end
                sum_Coulomb_k_ell_over_j = 0;
                
                partial_g_sum = 0;
                for m = 2 : 1 : M
                    B_t_m = B_t( :, m );
                    B_k_t_m = BB_component( B_t_m, k, N_e, dim );
                    BB_k_ell_t_m = B_k_t_m + ( 1 - ( m - 1 ) / M ) * x0_k + ( m - 1 ) / M * x0_ell;
                    partial_g_summand_1 = dt * ( BB_k_ell_t_m * B_k_t_m' ) / ( 2 * beta );
                    partial_Coulomb_inter_sum = 0;

                    partial_g_summand_2 = dt * partial_Coulomb_inter_sum;
                    partial_g_sum = partial_g_sum + partial_g_summand_1 + partial_g_summand_2;
                end
                a_k_ell = exp( -x0_diff_square / ( 2 * beta ) );
                V_tilde_sum = V_HO_sum + sum_Coulomb_k_ell_over_j;
                b_k_ell = exp( -V_tilde_sum );
                W_mat( k, ell ) = a_k_ell * b_k_ell;
                da_k_ell_dbeta = a_k_ell * x0_diff_square / ( 2 * beta^2 );
                dg_dbeta = partial_g_sum + V_tilde_sum / beta;
                db_k_ell_dbeta = b_k_ell * ( -dg_dbeta );
                W_deriv_mat( k, ell ) = da_k_ell_dbeta * b_k_ell + a_k_ell * db_k_ell_dbeta;
            end
        end
        
    end
        
        
        
        
        
    
    
    
    
    
    
    
    
    