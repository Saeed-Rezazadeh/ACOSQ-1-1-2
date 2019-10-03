function D_c = Channel_Distortion_step_3 (y_1 , y_2 , f , T , delta_u , Pr_z , codebook , b)
summation = 0 ;

for x_3 = 1 : 2
    for x_4 = 1 : 2
        x_3_4 = (x_3 - 1) * 2 + x_4 ;
        for y_3 = 1 : 2
            for y_4 = 1 : 2
                y_3_4 = (y_3 - 1) * 2 + y_4 ; 
                u_index = find (T (: , 5) == x_3_4) ;
                for u_i = 1 : length(u_index)
                    x_2 = T (u_index(u_i) , 2 + y_1) ;
                    summation = summation ...
                        + Pr_z(xor(b(x_2) - 1 , y_2 - 1) + 1 , xor(b(x_3) - 1 , y_3 - 1) + 1 ) ...
                        * Pr_z(xor(b(x_3) - 1 , y_3 - 1) + 1 , xor(b(x_4) - 1 , y_4 - 1) + 1 ) ...
                        * delta_u * f(u_index(u_i)) * (codebook (b(x_3_4)) - codebook(y_3_4)) ^ 2;
                end
            end
        end
    end
end
D_c = summation ;
end