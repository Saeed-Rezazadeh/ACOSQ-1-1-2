function D_c = Channel_Distortion_step_1 (f , T , numLevel , delta_u , Pr , codebook , b)
summation = 0 ; 

for i = 1 : numLevel
    for j = 1 : numLevel
        u_index = find (T (: , 2) == i) ;
        P_i = delta_u * sum (f(u_index)) ;
        summation = summation + Pr(j , b(i)) * P_i * (codebook (b(i)) - codebook(j)) ^ 2;
    end
end
D_c = summation ;
end