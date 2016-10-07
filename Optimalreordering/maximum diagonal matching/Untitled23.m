B = A30000;
reorder_file = fopen(sprintf('reorder max_min_dia.txt'),'wt+');

m = size(B, 2);
n = size(B, 1);
if m ~= n
    fprintf('Error! The input matrix is not square \n');
end

CUT = 5*max(norm(B, 1), norm(B, inf));
% C = B;
D = B;

Store_dia = zeros(1, n);

fixed_num = 1;
while fixed_num < m
%     fprintf(reorder_file, 'At step %d : \n', fixed_num);
    fprintf('At step %d : \n', fixed_num)
    % Find the cutoff in the diagonal with the eleminated matrix E
    if fixed_num == 1
        % Initial step
        E = B;
        [ cutoff, bisec_it ] = max_min_diagonal( E, 2*size(E, 1)+2, 1e-5, 1e-5, 50  );
    else
        % Steps after we change the diagonal elements of fixed rows
        [ cutoff, bisec_it ] = max_min_diagonal( E, 2*size(E, 1)+2, 1e-5, 1e-5, 50  );
    end
    % Reorder the matrix with new cutoff
    [ D, ~, Matching ] = max_weighted_match( D, cutoff );
    
    % Find the minimum row of the diagonal
    [dia_val, fix_row] = min(abs(diag(D)));
    fprintf('The smallest diagonal value is %f at row %d. \n', dia_val, fix_row);
    % Make sure possible swapping set exists.
    col = D(:, fix_row);
    if max(col) <= cutoff
        error('no element is larger than cutoff in %d column.\n', fix_row);
        break;
    end
    % Find all possible swaps with the row to be fixed
    correspond = zeros(2, n);
    correspond(1, :) = col';
    correspond(2, :) = D(fix_row, :);
    min_correspond = min(correspond);
    flag = 0;
    while flag == 0
        [max_value, max_index] = max(min_correspond);
        if max_value <= dia_val
            fprintf('No swap because the diagonal is already laregest. \n')
            break;
        end
        % Swap max_index row with the row we want to fix
        temp = D(max_index, :);
        D(max_index, :) = D(fix_row, :);
        D(fix_row, :) = temp;
        flag = 1;
        fprintf('Swap row %d and %d. \n', max_index, fix_row)
        % Check if the perfect matching exists if we make this swap
        % ???: should it always exist???
        CheckM = D;
        CheckM(fix_row, :) = [];
        CheckM(:, fix_row) = [];
    end
            % Update and store diagonal element
        dia_val = D(fix_row, fix_row);
        Store_dia(fix_row) = dia_val;
            % Penalte the diagonal element of the row we want to fix
            D(fix_row, fix_row) = -CUT;
            fprintf('The diagonal element at %d is %f. \n', fix_row, -CUT)
            % Resize our eliminated matrix E
            E(fix_row, :) = [];
            E(:, fix_row) = [];
%     [max_y, max_index_col] = max(correspond(2, :));
    fixed_num = fixed_num + 1;
end
C = D;
