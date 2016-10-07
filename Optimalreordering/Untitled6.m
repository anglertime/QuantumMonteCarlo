   

A = A12000; 

D = A;

max_each_col = max(D);
    [largest_dist, orbital] = min(max_each_col);
%     fprintf(reorder_file, 'For orbital %d \n', orbital);
    % Find the row number where min_value is located
