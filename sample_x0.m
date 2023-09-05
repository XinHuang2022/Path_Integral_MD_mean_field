    
    function [ x_0, p_x_density ] = sample_x0( beta, dim, N_e )
        
        d = N_e * dim;
        
		density_urv = rand;
		sigma_1 = 1 / sqrt( beta );
		sigma_2 = sqrt( beta );
		if( density_urv < 0.5 )
			x_0 = randn( d, 1 ) * sigma_1;
		else
			x_0 = randn( d, 1 ) * sigma_2;
		end
		p_x_density_1 = exp( -( x_0' * x_0 ) / ( 2 * sigma_1^2 ) ) / ( sqrt( 2 * pi ) * sigma_1 )^d;
		p_x_density_2 = exp( -( x_0' * x_0 ) / ( 2 * sigma_2^2 ) ) / ( sqrt( 2 * pi ) * sigma_2 )^d;
		p_x_density = ( p_x_density_1 + p_x_density_2 ) * 0.5;
        
    end
    
    
    
    