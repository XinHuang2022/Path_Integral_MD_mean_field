    
    function Z_value = Z_partition_a( n, beta, dim )
    
        if( n == 0 )
            Z_value = 1;
        elseif( n == 1 )
            Z_value = Z1_func( beta, dim );
        else
            Z_sum = 0;
            for i = 1 : 1 : n
                Z_sum = Z_sum + ( -1 )^( i + 1 ) * Z_partition_a( 1, i * beta, dim ) * Z_partition_a( n - i, beta, dim );    % for fermions
                % Z_sum = Z_sum + ( 1 )^( i + 1 ) * Z_partition( 1, i * beta, dim ) * Z_partition( n - i, beta, dim );    % for bosons
            end
            Z_value = Z_sum / n;
        end
    
    
    
    
    
    
    
    
    