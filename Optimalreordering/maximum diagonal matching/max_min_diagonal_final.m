function [ C, fixed_num ] = max_min_diagonal_final( B )
%Optimization. Maximize the smallest element in the diagonal of a square
%matrix by row permutations

% [Y,I] = min(X) returns the indices of the minimum values in vector I.
%     If the values along the first non-singleton dimension contain more
%     than one minimal element, the index of the first one is returned.

%REMARK:
%When three or even more particles are around an orbital and there are two
%or even more orbitals are almost alone, it is not good. The min max may
%not be easy to obtain.

maxminfile = fopen('reorder max_min_dia.txt','wt+');
m = size(B, 2);
n = size(B, 1);
if m ~= n
    fprintf('Error! The input matrix is not square \n');
end

% Initialization
CUT = 5*max(norm(B, 1), norm(B, inf));
D = B;
Store_dia = zeros(1, n);
fixed_num = 1;

%****************************************
%   Main Function
%****************************************

while fixed_num < m
    fprintf(maxminfile, '\n         At step %d : \n', fixed_num);
    fprintf('\n         At step %d : \n', fixed_num)
    
    %****************************************
    %   Step 1: Fix the first max-min diagonal
    %****************************************
    
    if fixed_num == 1
        % Initial step
        E = D;
        [ cutoff, bisec_it ] = max_min_diagonal( E, 2*size(E, 1)+2, 1e-5, 1e-5, 50  );
        
        % initial_cutoff will be used to check if we have obtained new cutoff in next step
        initialcutoff = cutoff;
        
        [ D, ~, Matching ] = max_weighted_match( D, cutoff ); % Reorder the matrix with cutoff
        % Find the minimum row of the diagonal
        [dia_val, fix_row] = min(abs(diag(D)));
        fprintf('The smallest diagonal value is %f at row %d. \n', dia_val, fix_row);
        fprintf(maxminfile, 'The smallest diagonal value is %f at row %d. \n', dia_val, fix_row);
        
        % This is to make sure possible swapping set exists.
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
                fprintf('NO swap necessary for row %d. \n', fix_row)
%                 fprintf(maxminfile, 'NO swap necessary for row %d. \n', fix_row);
                break;
            end
            % Swap max_index row with the row we want to fix
            temp = D(max_index, :);
            D(max_index, :) = D(fix_row, :);
            D(fix_row, :) = temp;
            flag = 1;
            fprintf('Swap row %d and %d. \n', max_index, fix_row)
            fprintf(maxminfile, 'Swap row %d and %d. \n\n\n', max_index, fix_row);
        end
        % Update and store diagonal element
        dia_val = D(fix_row, fix_row);
        Store_dia(fix_row) = dia_val;
        % Penalte the diagonal element of the row we want to fix
        D(fix_row, fix_row) = CUT;
%         fprintf('The diagonal element at %d is %f. \n', fix_row, CUT)
        % Resize our eliminated matrix E
        E(fix_row, :) = [];
        E(:, fix_row) = [];
        fixed_num = fixed_num + 1;
        rowsum = fix_row;
        
    elseif (fixed_num > 1)
        
        %****************************************
        %   Step 2: Keep fixing heuristically
        %****************************************
        
        % Find the new cutoff in the diagonal with the eleminated matrix E
        [ cutoff, bisec_it ] = max_min_diagonal( E, 2*size(E, 1)+2, 1e-5, 1e-5, 50  );
        flagout = 0;
        while flagout == 0
            
            if (cutoff < initialcutoff) %&& (dia_val > max(:, collllll))
                
                %****************************************
                %   Case 1: After we fixed the previous row, the new cutoff is not
                %   larger. No reordering needed.
                %****************************************
                % No reordering needed. Find the minimum and fix it!
                [dia_val, fix_row] = min(abs(diag(D)));
                fprintf('No reordering needed. Fix row %d with diagonal element %f. \n', fix_row, dia_val);
                fprintf(maxminfile, 'No reordering needed. Fix row %d with diagonal element %f. \n', fix_row, dia_val);
                rowsum = union(rowsum, fix_row); % Update the fixed rows set.
                Store_dia(fix_row) = dia_val;
                % Penalte the diagonal element of the row we want to fix
                D(fix_row, fix_row) = CUT;
                
                % Resize D to obtain the eliminated matrix E
                E = D;
                E(rowsum, :) = [];
                E(:, rowsum) = [];
                fixed_num = fixed_num + 1;
                break;
%             elseif (cutoff = initialcutoff)
                
            elseif (cutoff >= initialcutoff)
                
                %****************************************
                %   Case 2: Reordering needed for new and larger cutoff
                %****************************************
                
                % Reorder the matrix with new cutoff
                [ D, ~, Matching ] = max_weighted_match( D, cutoff );
                initialcutoff = cutoff;  % Update initialcutoff
                
                % Find the minimum row of the diagonal
                [dia_val, fix_row] = min(abs(diag(D)));
%                 fprintf('The smallest diagonal value is %f at row %d. \n', dia_val, fix_row);
%                 fprintf(maxminfile, 'The smallest diagonal value is %f at row %d. \n', dia_val, fix_row);
                % This is to make sure possible swapping set exists.
                col = D(:, fix_row);
                if max(col) < cutoff
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
%                         fprintf('NO swap necessary for row %d. \n\n', fix_row)
%                         fprintf(maxminfile, 'NO swap necessary for row %d. \n\n', fix_row);
                        break;
                    end
                    % Swap max_index row with the row we want to fix
                    temp = D(max_index, :);
                    D(max_index, :) = D(fix_row, :);
                    D(fix_row, :) = temp;
                    flag = 1;
                    fprintf('Swap row %d and %d. \n\n\n', max_index, fix_row)
                end
                % Update and store diagonal element
                dia_val = D(fix_row, fix_row);
                Store_dia(fix_row) = dia_val;
                % Penalte the diagonal element of the row we want to fix
                D(fix_row, fix_row) = CUT;
%                 fprintf('Row %d is to be fixed; diagonal ele is %f. \n', fix_row, dia_val)
%                 fprintf(maxminfile, 'Row %d is to be fixed; diagonal ele is %f. \n', fix_row, dia_val);
                rowsum = union(rowsum, fix_row); % Update the fixed rows set.
                % Resize D to obtain the eliminated matrix E
                E = D;
                E(rowsum, :) = [];
                E(:, rowsum) = [];
                
                fixed_num = fixed_num + 1;
                flagout = 1;
            end
        end
        
    end
end
C = D;
first = find(Store_dia == 0);
Store_dia(first) = D(first, first);
C(diag(true(m,1))) = Store_dia;
fclose(maxminfile);