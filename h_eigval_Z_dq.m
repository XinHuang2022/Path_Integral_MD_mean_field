    
    function h_dq = h_eigval_Z_dq( N_e, beta, dim )
        
        Z_beta_1 = Z_partition_a( N_e, beta, dim );
        threshold = 10^(-7);
        delta_beta = 0.1;
        dq_diffence = 1;
        beta_2 = beta + delta_beta;
        Z_beta_2 = Z_partition_a( N_e, beta_2, dim );
        h_a = -( log( Z_beta_2 ) - log( Z_beta_1 ) ) / ( delta_beta );
        iter_count = 1;
        while( dq_diffence > threshold )
            delta_beta = delta_beta / 2;
            beta_2 = beta + delta_beta;
            Z_beta_2 = Z_partition_a( N_e, beta_2, dim );
            h_b = -( log( Z_beta_2 ) - log( Z_beta_1 ) ) / ( delta_beta );
            dq_diffence = abs( h_a - h_b ) / abs( h_b );
            h_a = h_b;
            iter_count = iter_count + 1;
        end
        h_dq = h_b;
        % iter_count
    
    end
    
    
    
    
    
    
    
    