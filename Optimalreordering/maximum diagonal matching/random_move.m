function [movingP, move] = random_move(N, sigma, dim)

movingP = ceil(rand(1)*N);
move = randn(1,dim)*sigma;

