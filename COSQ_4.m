function [SDR , Distortion , T , codebook] = COSQ_4(f  , y_1 , y_2 , Pr_z , T , codebook , delta )
y_1_2 = (y_1 - 1) * 2 + y_2 ;
FileID = fopen ('Results.txt' , 'a') ;

D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    T_u = zeros(length(T) , 1) ; 
    parfor u_index = 1 : length(T) 
        d_4 = zeros(4 , 1) ; 
        summation = 0 ; 
        
        u = T(u_index , 1) ; 
        x_2 = T(u_index , 2 + y_1) ; 
        for x_3 = 1 : 2 
            for x_4 = 1 : 2 
                x_prime = (x_3 - 1) * 2 + x_4 ; 
                for y_3 = 1 : 2 
                    for y_4 = 1 : 2 
                        y_prime = (y_3 - 1) * 2 + y_4 ; 
                        summation = summation + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) * ... 
                            Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) * (u - codebook(y_prime)) ^ 2 ; 
                    end 
                end 
                d_4(x_prime) = summation ; 
                summation = 0 ; 
            end 
        end
        [~ , partition_index] = min(d_4) ; 
        T_u(u_index , 1) = partition_index ; 
    end 
    T(: , 5) = T_u ; 
    %% Optimal Centroids
    for y_3 = 1 : 2 
        for y_4 = 1 : 2 
            y_prime = (y_3 - 1) * 2 + y_4 ; 
            numerator = 0 ; 
            denominator = 0 ; 
            for x_3 = 1 : 2 
                for x_4 = 1 : 2 
                    x_prime = (x_3 - 1) * 2 + x_4 ; 
                    u_index = find(T(: , 5) == x_prime ) ; 
                    
                    for u_i = 1 : length(u_index) 
                        x_2 = T(u_index(u_i) , 2 + y_1) ; 
                        
                        numerator = numerator + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) ...
                            * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                            * T(u_index(u_i) , 1) * f(u_index(u_i)) ;
                        
                        denominator = denominator + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) ...
                            * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                            * f(u_index(u_i)) ;
                    end 
                end 
            end
            codebook(y_prime) = numerator / denominator ; 
        end 
    end 
    %% Distortion
    D(2) = distortion_4(f , y_1 , y_2 , codebook , delta , Pr_z , T) ;
    fprintf (FileID , 'Overall D_4 = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_4 = %f\n' , SDR) ;
fclose (FileID) ;
end