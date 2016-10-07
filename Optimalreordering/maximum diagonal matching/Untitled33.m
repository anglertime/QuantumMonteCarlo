% A = dlmread('12000Matrix.txt');
% A = A_shifted; 
% A = [4,3,2,1;
%      1,4,3,2;
%      2,1,4,3;
%      2.5,3.5,4.5,1.5];
% [Matching, mat_minimum_weight, Cost] = Hungarian(-A);
% A_maximum_weight = -mat_minimum_weight;
% zero = zeros(4, 4);
% % G = [zero A; A' zero];
% % [ indices, cost ] = min_perfect_matching( G );
% cutoff = 3;
% % [i,j,s] = find(abs(A>cutoff).*A);
% % A_cut = sparse(i,j,s);
% % A_cut = full(A_cut);
% % A_cut = abs(A>cutoff).*A;
% A_cut = abs(A>cutoff);
% G = [zero A_cut; A_cut' zero];
% [ indices, ~, match ] = min_perfect_matching( G )
% [A_do] = make_diagonal_dominant(A);

% cutoff = 8e-3;
% A_adjacent = abs(A>cutoff);
% n = size(A_adjacent);
% % zero_A = zeros(size(A_adjacent));
% m = 2*n + 2;
% A_construct = zeros(m);
% for i = 2:n+1
%     A_construct(1, i) = 1;
% end
% for j = 1:n
%     A_construct( 2*n + 2 - n:2*n + 2 - 1,2*n + 2 ) = 1;
% end
% A_construct( 2:2+n-1, n+2:n+2 + (n-1)) = A_adjacent;
% % A_construct( n+2:n+2 +(n-1), 2:2+n-1) = A_adjacent';
% % [i,j,s] = find(A_construct>0);
% % A_construct = sparse(i,j,s,m-1,m-1, 10*n);
% A_construct = sparse(A_construct);
% 
% [flowval cut R F] = max_flow(A_construct, 1, 82);
% matching = F( 2:2+n-1, n+2:n+2 + (n-1));
% [part, orb] = find(matching>0);
% 
% A_opt = A;
% for i = 1:n
%     % go through the matching
%     % For orbital orb(i), swap part(i) into the diagonal.
%     
% %     for innerCount = 1:n
% %     end
%     temp = A_opt(orb(i), :);
%     A_opt(orb(i), :) = A_opt(part(i), :);
%     A_opt(part(i), :) = temp;
%     
%     % Swap particles orb(i) and part(i)
%     temp = part(part == orb(i));
%     part(part == orb(i)) = part(i);
% end
% A = dlmread('12000Matrix.txt');
% A = A_shifted; 
B = A;
tic;
[ cutoff, bisec_it ] = max_min_diagonal( B, 2*size(B, 1)+2, 1e-6 );
% particle = dlmread('30000particle.txt');

% [ B_opt, cost, newParticle, Matching ] = max_weighted_match( A, cutoff );
t = toc



% [ C_do ] = make_diagonal_dominant( A );
% [ C, swap_nr ] = max_min_diagonal_final( A );
% cond(B_opt_maxweit)
% cond(A)
% cond(A_opt)
% plot(real(eig(A_opt)),imag(eig(A_opt)), '*k')
% hold on
% plot(real(eig(A)),imag(eig(A)), '*b')
% grid on
% hold on
% plot(real(eig(B_opt)),imag(eig(B_opt)), '*m')
% hold on
% plot(real(eig(C_do)),imag(eig(C_do)), '*r')