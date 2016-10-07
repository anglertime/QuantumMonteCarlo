function [ C, swap_nr ] = mix_max_min_diagonal_final( B, tol )
%Optimization. Maximize the smallest element in the diagonal of a square
%matrix by row permutations

% [Y,I] = min(X) returns the indices of the minimum values in vector I.
%     If the values along the first non-singleton dimension contain more
%     than one minimal element, the index of the first one is returned.

%REMARK:
%When three or even more particles are around an orbital and there are two
%or even more orbitals are almost alone, it is not good. The min max may
%not be easy to obtain.

reorder_file = fopen(sprintf('reorder max_min_dia.txt'),'wt+');

n1 = size(B, 2);
n = size(B, 1);
if n1 ~= n
    fprintf('Error! The input matrix is not square \n');
end

CUT = max(max(B)) + 10;
C = B;
D = B;
k = 1;
particles = [];
orbitals = [];

%%
while k < n-1
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
    max_it = 50;
    kk = 0; start = min_val; terminal = max_val; cutoff_m = (start + terminal)/2;
    while (kk < max_it)       %     fprintf('\n\nBisection step %d:  \n[A, B] = [%2.6f, %2.6f], cutoff_m = %2.6f \n', k + 1, start, terminal, cutoff_m);
        B_adjacent = abs( B>cutoff_m );
        B_construct( 2:2+n-1, n+2:n+2 + (n-1)) = B_adjacent;
        B_construct = sparse(B_construct);
        [max_flow_value, ~, ~, F] = max_flow(B_construct, 1, size_nrm);
        if max_flow_value == n
            start = cutoff_m;
        elseif max_flow_value < n
            terminal = cutoff_m;
        end
        cutoff_m = (start + terminal)/2;        %     fprintf('[A, B] = [%2.6f, %2.6f], cutoff_m = %2.6f. max_flow is %d \n', start, terminal, cutoff_m, max_flow_value);
        if (cutoff_m > 1e-8) && (abs(start - terminal) < tol )
            break;
        end
        kk = kk + 1;
    end
    cutoff_m = start;
    B_adjacent = abs( B>cutoff_m );
    B_construct( 2:2+n-1, n+2:n+2 + (n-1)) = B_adjacent;
    B_construct = sparse(B_construct);
    [max_flow_value, ~, ~, F] = max_flow(B_construct, 1, size_nrm);
    
    fprintf('\nBisection steps %d, cutoff value %2.12f, |E_m| is %d \n\n', kk, cutoff_m, max_flow_value);
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
        temp = D(orb(i), :);
        D(orb(i), :) = D(part(i), :);
        D(part(i), :) = temp;
        % Swap particles orb(i) and part(i)
        %     temp = part(part == orb(i));
        part(part == orb(i)) = part(i);
    end
    
    k = k + 1;
end

end

