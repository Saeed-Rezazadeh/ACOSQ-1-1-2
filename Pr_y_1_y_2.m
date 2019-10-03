function Probability_y_1_2 = Pr_y_1_y_2( y_1, y_2 , Pr_1 , Pr_z , f  , T , delta)
summation = 0 ;
for x_1 = 1 : 2
    for x_2 = 1 : 2
        u_index = find (T(: , 2) == x_1 & T(: , 2 + y_1) == x_2) ;
        for u_i = 1 : length(u_index)
            summation = summation + delta * Pr_1(x_1 , y_1) ...
                *  Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) * f(u_index(u_i)) ;
        end
    end
end
Probability_y_1_2 = summation ;
end