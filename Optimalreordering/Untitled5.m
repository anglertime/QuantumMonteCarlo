tic
% [Matching, Cost, num_y, num_x, x_con, y_con, P_cond, P_size ] = Hungariantest(-A);
[assignment,cost] = munkres(-A);
t = toc