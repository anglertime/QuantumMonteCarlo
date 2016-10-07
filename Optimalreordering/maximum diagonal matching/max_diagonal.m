function [ B ] = max_diagonal( B, cutoff_tol )
%Optimization. Maximize the smallest element in the diagonal of a square
%matrix by row permutations

% [Y,I] = min(X) returns the indices of the minimum values in vector I.
%     If the values along the first non-singleton dimension contain more
%     than one minimal element, the index of the first one is returned.

%REMARK:
%When three or even more particles are around an orbital and there are two
%or even more orbitals are almost alone, it is not good. The min max may
%not be easy to obtain.

file = fopen(sprintf('max_min_dia.txt'),'wt+');
n1 = size(B, 2);
n = size(B, 1);
if n1 ~= n
    fprintf('Error! The input matrix is not square \n');
end

C = B;
pointer = zeros(1, n1);
for i = 1:n1
    pointer(i) = i;
end
% delete = zeros(1, n1);
% store = zeros(n1, 2*n1);
for i = 1:n1-1
    % Find threshold for C
    [ cutoff, k ] = bisect_threshold( C, 50, 2*size(C, 1)+2, cutoff_tol );
    % Find the smallest diagonal element
    [dmin, col_index] = min(diag(C));
    % Fix it with the maximal diagonal choice in column(row_index)
    row = C(col_index, :);
    col = C(:, col_index);
    fesible = [row; col'];
    fesible = min(fesible);
    for k = 1:length(row)
        if k == col_index
            fesible(k) = 0;
        else
            if fesible(k) > cutoff
                fesible(k) = 0;
            end
        end
    end
    [~, cut_index] = max(fesible);
    % Swap row(cut_index) and row(col_index)
    temp = C(cut_index, :);
    C(cut_index, :) = C(col_index, :);
    C(col_index, :) = temp;
    % Swap in B
    B_row2 = find(pointer == cut_index);
    B_row1 = find(pointer == col_index);
    temp = B(B_row2, :);
    B(B_row2, :) = B(B_row1, :);
    B(B_row1, :) = temp;
%     % Store row and columns
%     store(:, 2*i) = C(:, col_index);
%     store(:, 2*i-1) = C(col_index, :)';
%     % Store row number
%     index = find(pointer == col_index);
%     delete(i) = index;
    % Resize C
    C(:, col_index) = []; C(col_index, :) = []; 
    % New matching for row/column
    for j = B_row1:n1
        if j == B_row1
            pointer(j) = 0;
        else
            pointer(j) = pointer(j) - 1;
        end
    end
    fprintf(file, 'At step %d, Cutoff is %f. We swap row %d and row %d. \n', i, cutoff, B_row2, B_row1);
    fprintf(file, 'We fix row %d. \n', B_row1);
end
%%
fclose(file);
end

