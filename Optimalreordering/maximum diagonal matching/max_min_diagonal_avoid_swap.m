function [ C, swap_nr ] = max_min_diagonal_avoid_swap( B )
%Optimization. Maximize the smallest element in the diagonal of a square
%matrix by row permutations

% [Y,I] = min(X) returns the indices of the minimum values in vector I.
%     If the values along the first non-singleton dimension contain more
%     than one minimal element, the index of the first one is returned.
reorder_file = fopen(sprintf('reorder max_min_dia.txt'),'wt+');

m = size(B, 2);
n = size(B, 1);
if m ~= n
    fprintf('Error! The input matrix is not square \n');
end

flag = 0;
C = B;
D = B;
k = 0;
swap_row = [];

% while (flag == 0) && (k < 5)
while (flag == 0)
%     fprintf('At step %d \n', k);
    [min_dia, row_index] = min(diag(D));
    fprintf('The minimum is %2.4f in row %d \n', min_dia, row_index);
    row = D(row_index, :);
    col = D(:, row_index);
    % col should not be zeros because there must be some particles around
    % that orbital, otherwise it would be a zero column then the matrix is
    % absolutely singular.
%     if max(row) <= min_dia && max(col) <= min_dia
%         flag = 1;  
%         break;
%     end
    if max(col) <= min_dia
        flag = 1;  
        break;
    end
    correspond = zeros(2, n);
    for i = 1:n
        correspond(:, i) = [row(i); col(i)];
    end
    correspond(:, row_index) = [-1, -1]';
    min_correspond = min(correspond);
    [max_value, max_index] = max(min_correspond);
    if max_value < min_dia
        flag = 1;
        break;
    end
    for i = 1:n
        if min(correspond(:, i)) < min_dia
            correspond(:, i) = [-1; -1];
        end
    end
    
    [max_y, max_index_col] = max(correspond(2, :));
    
    % Swap row_index and max_index_col
    fprintf(reorder_file, '\n Swap row %d and %d \n', row_index, max_index_col);
    temp = C(max_index_col, :);
    C(max_index_col, :) = C(row_index, :);
    C(row_index, :) = temp;
    k = k + 1;
    D = C;
    
    temp_dia = D(row_index, row_index);
    D(row_index, :) = zeros(1, n);
    D(row_index, row_index) = temp_dia;


    [min_dia, row_index] = min(diag(D));
    fprintf('The minimum is %2.4f in row %d \n', min_dia, row_index);
    row = D(row_index, :);
    col = D(:, row_index);
    % col should not be zeros because there must be some particles around
    % that orbital, otherwise it would be a zero column then the matrix is
    % absolutely singular.
%     if max(row) <= min_dia && max(col) <= min_dia
%         flag = 1;  
%         break;
%     end
    if max(col) <= min_dia
        flag = 1;  
        break;
    end
    correspond = zeros(2, n);
    for i = 1:n
        correspond(:, i) = [row(i); col(i)];
    end
    correspond(:, row_index) = [-1, -1]';
    min_correspond = min(correspond);
    [max_value, max_index] = max(min_correspond);
    if max_value < min_dia
        flag = 1;
        break;
    end
    for i = 1:n
        if min(correspond(:, i)) < min_dia
            correspond(:, i) = [-1; -1];
        end
    end
    
    [max_y, max_index_col] = max(correspond(2, :));
    
    % Swap row_index and max_index_col
    fprintf(reorder_file, '\n Swap row %d and %d \n', row_index, max_index_col);
    temp = C(max_index_col, :);
    C(max_index_col, :) = C(row_index, :);
    C(row_index, :) = temp;
    k = k + 1;
    D = C;
    % avoid swap A and B cyclicaly
%     swap_row = union(swap_row, row_index);
%     temp_dia = D(row_index, row_index);
%     D(row_index, :) = zeros(1, n);
%     D(row_index, row_index) = temp_dia;
end
swap_nr = k;
fclose(reorder_file);
end

