function Distortion = distortion_2(f ,  y_1 , codebook , delta , Pr_z , T)
summation = 0 ;
parfor x_2 = 1 : 2
    for y_2 = 1 : 2
        u_index = find (T(: , 3) == x_2) ;
        for u_i = 1 : length(u_index)
            x_1 = T(u_index(u_i) , 2) ; 
            summation = summation + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 )...
                * delta * f(u_index(u_i)).* (T(u_index(u_i) , 1) - codebook(y_2)) .^ 2 ;
        end
    end
end
Distortion = summation ;
end