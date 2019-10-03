function D_c = Channel_Distortion_step_2 (y_1 , f , T , delta_u , Pr_z , codebook , b)
summation = 0 ;

for x_2 = 1 : 2
    for y_2 = 1 : 2
        u_index = find (T (: , 3) == x_2) ;
        for u_i = 1 : length(u_index)
            x_1 = T (u_index(u_i) , 2) ; 
            summation = summation ...
                + Pr_z(xor(b(x_1) - 1 , y_1 - 1) + 1 , xor(b(x_2) - 1 , y_2 - 1) + 1 ) ...
                * delta_u * f(u_index(u_i)) * (codebook (b(x_2)) - codebook(y_2)) ^ 2;
        end
    end
end
D_c = summation ;
end