function [ B_opt, cutoff_m, k ] = bisect_threshold( B, max_it, size_nrm, tol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% size_nrm must be an integer
% k: number of iterations of bisection method
%%
% Initialize

n1 = size(B, 2);
n = size(B, 1);
if n1 ~= n
    fprintf('Error! The input matrix is not square \n');
end
% Initialize for graph adjacent matrix
m = 2*n + 2;
B_construct = zeros(m);
for i = 2:n+1
    B_construct(1, i) = 1;
end
for j = 1:n
    B_construct( 2*n + 2 - n:2*n + 2 - 1,2*n + 2 ) = 1;
end

max_val = max(max(B));
min_val = 0.0;
k = 0; start = min_val; terminal = max_val; cutoff_m = (start + terminal)/2;

while (k < max_it)    
    B_adjacent = abs( B>cutoff_m );
    B_construct( 2:2+n-1, n+2:n+2 + (n-1)) = B_adjacent;
    B_construct = sparse(B_construct);
    [max_flow_value, ~, ~, ~] = max_flow(B_construct, 1, size_nrm);
    
    if max_flow_value == n
        start = cutoff_m;
    elseif max_flow_value < n
        terminal = cutoff_m;
    end
    cutoff_m = (start + terminal)/2;
    if (cutoff_m > 1e-8) && (abs(start - terminal) < tol )
        break;
    end
    k = k + 1;
end
cutoff_m = start;
B_adjacent = abs( B>cutoff_m );
B_construct( 2:2+n-1, n+2:n+2 + (n-1)) = B_adjacent;
B_construct = sparse(B_construct);
[max_flow_value, ~, ~, F] = max_flow(B_construct, 1, size_nrm);
fprintf('\nBisection steps %d, cutoff value %2.12f, |E_m| is %d \n\n', k, cutoff_m, max_flow_value);
%%
matching = F( 2:2+n-1, n+2:n+2 + (n-1));
[part, orb] = find( matching > cutoff_m );
B_opt = B;
for i = 1:n
    % go through the matching
    % For orbital orb(i), swap part(i) into the diagonal. 
%     for innerCount = 1:n
%     end
    temp = B_opt(orb(i), :);
    B_opt(orb(i), :) = B_opt(part(i), :);
    B_opt(part(i), :) = temp;
    % Swap particles orb(i) and part(i)
%     temp = part(part == orb(i));
    part(part == orb(i)) = part(i);
end
end

