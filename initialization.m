function [T , indexed_T , width] = initialization (sigma , u , alpha , resolution)
% Generate the training set based on the source parameters 
T = sigma * randn (alpha , 1) + u;
MAX = 4 * sigma ; 
MIN = -4 * sigma ; 
width = (MAX - MIN) / resolution ;
indexed_T = zeros (resolution + 1, 2) ;
indexed_T(: , 1) = MIN : width : MAX ;
end
