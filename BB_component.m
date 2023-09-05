    
    function BB_k = BB_component( BB_vec, k, N_e, dim )
        
        if( dim == 3 )
            BB_k = [ BB_vec( k, 1 ), BB_vec( N_e + k, 1 ), BB_vec( 2 * N_e + k, 1 ) ];
        elseif( dim == 2 )
            BB_k = [ BB_vec( k, 1 ), BB_vec( N_e + k, 1 ) ];
        elseif( dim == 1 )
            BB_k = BB_vec( k, 1 );
        end
        
    end
        
        
    
    
    
    