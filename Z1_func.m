    
    function Z1 = Z1_func( beta, dim )
        
        %{
        N_eigval = 10000;
		sum_1 = 0;
		for j_a = 0 : 1 : N_eigval
            E_j = ( j_a + 0.5 );
            sum_1 = sum_1 + exp( -beta * E_j );
        end
		Z1 = sum_1^dim;
		%}
        
        Z1 = ( exp( -beta / 2 ) / ( 1 - exp( -beta ) ) )^dim;
        
	end
    
    %{
    for j_a = 1 : 1 : N_eigval
            for j_b = 1 : 1 : N_eigval
                for j_c = 1 : 1 : N_eigval
                    E_j = 
	    end
    %}
    
    
    
    
    
    
    