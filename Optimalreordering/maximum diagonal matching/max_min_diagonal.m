function [ cutoff, bisec_it] = max_min_diagonal( B, size_BGraph, cutoff_tol, bi_tol, max_it  )
% Maximize the smallest element in the diagonal of B. cutoff is the max-min solution
%
% Use max flow method to determine if the perfect matching exists for a
% matrix. Employ bisection method to obtain the optimal value.
%
% Cutoff: Round cutoff to cutoff_tol digits to get a bit more choices of matching.
% size_B must be an integer equals 2n+2, where n^2 is the dimension of B.
% bisec_it: number of iterations of bisection method
% bi_tol : tol in the bisection process.
% cutoff_tol : tol we want to round the cutoff value. It should be 1e-2,
% 1e-3, 1e-4, etc.
%%
% Initialize
n1 = size(B, 2);
n = size(B, 1);
if n1 ~= n
    fprintf('Error! The input matrix is not square \n');
end
% Initialize adjacent matrix for directed graph
m = 2*n + 2;
B_construct = zeros(m);
B_construct(1, 2:n+1) = 1;
B_construct( 2*n+2 - n:2*n+2 -1, 2*n+2 ) = 1;
%%
max_val = max(max(B));
min_val = 0.0;
% max_it = 50;
bisec_it = 1; start = min_val; terminal = max_val; cutoff = (start + terminal)/2;
% Bisection algorithm
while (bisec_it < max_it) && (abs(start - terminal) > bi_tol )
    B_adjacent = abs( B>cutoff );
    B_construct( 1+ 1:1+ n, n+1+ 1:n+1+ n) = B_adjacent;
    B_construct = sparse(B_construct);
    
    % Push_relabel is the default algorithm. Edmunds-Karp and Kolmogorov
    % algorithms are also available.
    [max_flow_value, ~, ~, ~] = max_flow(B_construct, 1, size_BGraph);
%     [max_flow_value, ~, ~] = graphmaxflow(B_construct, 1, size_B); 
    
    if max_flow_value == n % Perfect matching exists
        start = cutoff;
    elseif max_flow_value < n % Perfect matching does not exist
        terminal = cutoff;
    elseif max_flow_value > n
        fprintf('Error! The max flow is more than the number of particles. \n');
        error('The max flow is more than the number of particles.');
    end
    cutoff = (start + terminal)/2;
%     if (cutoff < 1e-8) && (abs(start - terminal) < bi_tol )
%         break;
%     end
    bisec_it = bisec_it + 1;
end
cutoff = start; % Set cutoff equal to the start because the perfect matching may not exist for (start + terminal)/2;
cutoff = fix(cutoff/cutoff_tol )*cutoff_tol ; % Round cutoff to cutoff_tol  digits to get a bit more choices of matching.
B_adjacent = abs( B>cutoff );
B_construct( 2:2+n-1, n+2:n+2 + (n-1)) = B_adjacent;
B_construct = sparse(B_construct);

[max_flow_value, ~, ~, FlowMatrix] = max_flow(B_construct, 1, size_BGraph);
% [max_flow_value, FlowMatrix, ~] = graphmaxflow(B_construct, 1, size_B); 

fprintf('Bisection steps %d, cutoff value %2.12f, max flow value is %d \n', bisec_it, cutoff, max_flow_value);
%%

end