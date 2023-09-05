    
    function [ W_mat, W_deriv_mat ] = W_matrices_mod_V2( x_0, B_t, beta, dt, N_e, M, dim, lambda )
    
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
                for j = 1 : 1 : N_e
                    if( ( j ~= k ) )
                        x0_j = x0_vectors( j, : );
                        if( j ~= ell )
                            nu_j = j;
                        else
                            nu_j = k;
                        end

                        x0_nu_j = x0_vectors( nu_j, : );
                        V_Coulomb_integ_k_ell_j = ( lambda / ( 2 * norm( x0_k - x0_j ) ) + lambda / ( 2 * norm( x0_ell - x0_nu_j ) ) ) * dt / 2;
                        for m = 2 : 1 : M
                            B_t_m = B_t( :, m );
                            B_k_t_m = BB_component( B_t_m, k, N_e, dim );
                            B_j_t_m = BB_component( B_t_m, j, N_e, dim );
                            x_k_ell_t_m = B_k_t_m + ( 1 - ( m - 1 ) / M ) * x0_k + ( ( m - 1 ) / M ) * x0_ell;
                            x_j_nu_j_t_m = B_j_t_m + ( 1 - ( m - 1 ) / M ) * x0_j + ( ( m - 1 ) / M ) * x0_nu_j;
                            distance_k_ell_j_nu_j = norm( x_k_ell_t_m - x_j_nu_j_t_m );
                            V_Coulomb_k_ell_j_nu_j_summand = lambda / ( 2 * distance_k_ell_j_nu_j );
                            V_Coulomb_integ_k_ell_j = V_Coulomb_integ_k_ell_j + V_Coulomb_k_ell_j_nu_j_summand * dt;
                        end
                        sum_Coulomb_k_ell_over_j = sum_Coulomb_k_ell_over_j + V_Coulomb_integ_k_ell_j;
                    end
                end
                
                partial_g_sum = 0;
                for m = 2 : 1 : M
                    B_t_m = B_t( :, m );
                    B_k_t_m = BB_component( B_t_m, k, N_e, dim );
                    BB_k_ell_t_m = B_k_t_m + ( 1 - ( m - 1 ) / M ) * x0_k + ( m - 1 ) / M * x0_ell;
                    partial_g_summand_1 = dt * ( BB_k_ell_t_m * B_k_t_m' ) / ( 2 * beta );
                    partial_Coulomb_inter_sum = 0;
                    
                    for j = 1 : 1 : N_e
                        if j == k
                            continue
                        else
                            x0_j = x0_vectors( j, : );
                            B_j_t_m = BB_component( B_t_m, j, N_e, dim );
                            if( j == ell )
                                nu_j = k;
                            else
                                nu_j = j;
                            end
                            x0_nu_j = x0_vectors( nu_j, : );
                            B_j_nu_j_t_m = B_j_t_m + ( 1 - ( m - 1 ) / M ) * x0_j + ( m - 1 ) / M * x0_nu_j;
                            BB_diff_t_m = BB_k_ell_t_m - B_j_nu_j_t_m;
                            gradient_BB_diff = ( B_k_t_m - B_j_t_m ) / ( 2 * beta );
                            partial_Coulomb_inter_summand = -( BB_diff_t_m * gradient_BB_diff' ) * lambda / ( 2 * norm( BB_diff_t_m )^3 );
                            partial_Coulomb_inter_sum = partial_Coulomb_inter_sum + partial_Coulomb_inter_summand;
                        end
                    end
                    
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
        
        
        
        
        
    
    
    
    
    
    
    
    
    