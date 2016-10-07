
global sideln constK

datafile = fopen(sprintf('inexact -OutputK%d.txt',constK),'wt+');
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
level_1 = 0.5001*sideln;
ele_level_1 = (1/3.5449024699805363) * exp(-constK*level_1^2);
level_2 = 1*sideln;
ele_level_2 = (1/3.5449024699805363) * exp(-constK*level_2^2);
B = A;
% N = size(B, 1);
par = 1:13;
i = 1;
swap = [];


    fprintf(datafile, 'At step %d; \n', par(1));
    [row, local_orb] = find(B(par(1), :) > ele_level_1);
    local_orb = setdiff(local_orb, swap);
    if length(local_orb) == 1
        % Swap the value in B matrix for those two columns
        if local_orb == par(1)
            fprintf(datafile, 'No swap. Same: %d and %d \n', par(1), local_orb);
            swap = union(swap, par(1));
        else
            for innerCount = 1:N
                innerValue = B(innerCount, par(1));
                B(innerCount, par(1)) = B(innerCount, local_orb);
                B(innerCount, local_orb) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', par(1), local_orb);
            swap = union(swap, par(1));
%             swap = union(swap, local_orb);
        end
        par = setdiff(par, par(1));
        
    elseif length(local_orb) == 2
        fprintf(datafile, 'There are two orbitals \n');
        [max_a, indexa] = max(abs(B(1:N, local_orb(1))));
        [max_b, indexb] = max(abs(B(1:N, local_orb(2))));
        if max_a < max_b
            for innerCount = 1:N
                innerValue = B(innerCount, par(1));
                B(innerCount, par(1)) = B(innerCount, local_orb(1));
                B(innerCount, local_orb(1)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', par(1), local_orb(1));
            swap = union(swap, par(1));
%             swap = union(swap, local_orb(1));
            par = setdiff(par, par(1));
            % Assign the other particle indexb to local_orb(2)
            for innerCount = 1:N
                innerValue = B(innerCount, indexb);
                B(innerCount, indexb) = B(innerCount, local_orb(2));
                B(innerCount, local_orb(2)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', indexb, local_orb(2));
            swap = union(swap, indexb);
%             swap = union(swap, local_orb(2));
            par = setdiff(par, indexb);
        else
            for innerCount = 1:N
                innerValue = B(innerCount, par(1));
                B(innerCount, par(1)) = B(innerCount, local_orb(2));
                B(innerCount, local_orb(2)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', par(1), local_orb(2));
            swap = union(swap, par(1));
%             swap = union(swap, local_orb(2));
            par = setdiff(par, par(1));
            % Assign the other particle indexa to local_orb(1)
            for innerCount = 1:N
                innerValue = B(innerCount, indexa);
                B(innerCount, indexa) = B(innerCount, local_orb(1));
                B(innerCount, local_orb(1)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', indexb, local_orb(1));
            swap = union(swap, indexa);
%             swap = union(swap, local_orb(1));
            par = setdiff(par, indexa);
        end
    elseif length(local_orb) > 2
        fprintf(datafile, 'Error! There are more than two orbitals inside the level 1 of the particle \n');
        break;
    elseif isempty(local_orb) == 1
        fprintf(datafile, 'Zero orbital inside the level 1 of the particle %d after cutoof \n  Then we use level 2 cutoff. \n \n', par(1));
        C = B;
        for inner = 1:length(swap)
            C(par(1), swap(inner)) = 0.0;
        end
        [maxValue, indexc] = max(abs(C(par(1), :)));
        % Assign indexc to par(1)
        for innerCount = 1:N
            innerValue = B(innerCount, indexc);
            B(innerCount, indexc) = B(innerCount, par(1));
            B(innerCount, par(1)) = innerValue;
        end
        fprintf(datafile, 'So, swap %d and %d \n', indexc, par(1));
        swap = union(swap, indexc);
        par = setdiff(par, par(1));
    end         
    
    
    
    fprintf(datafile, 'At step %d; \n', par(1));
    [row, local_orb] = find(B(par(1), :) > ele_level_1);
    local_orb = setdiff(local_orb, swap);
    if length(local_orb) == 1
        % Swap the value in B matrix for those two columns
        if local_orb == par(1)
            fprintf(datafile, 'No swap. Same: %d and %d \n', par(1), local_orb);
            swap = union(swap, par(1));
        else
            for innerCount = 1:N
                innerValue = B(innerCount, par(1));
                B(innerCount, par(1)) = B(innerCount, local_orb);
                B(innerCount, local_orb) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', par(1), local_orb);
            swap = union(swap, par(1));
%             swap = union(swap, local_orb);
        end
        par = setdiff(par, par(1));
        
    elseif length(local_orb) == 2
        fprintf(datafile, 'There are two orbitals \n');
        [max_a, indexa] = max(abs(B(1:N, local_orb(1))));
        [max_b, indexb] = max(abs(B(1:N, local_orb(2))));
        if max_a < max_b
            for innerCount = 1:N
                innerValue = B(innerCount, par(1));
                B(innerCount, par(1)) = B(innerCount, local_orb(1));
                B(innerCount, local_orb(1)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', par(1), local_orb(1));
            swap = union(swap, par(1));
%             swap = union(swap, local_orb(1));
            par = setdiff(par, par(1));
            % Assign the other particle indexb to local_orb(2)
            for innerCount = 1:N
                innerValue = B(innerCount, indexb);
                B(innerCount, indexb) = B(innerCount, local_orb(2));
                B(innerCount, local_orb(2)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', indexb, local_orb(2));
            swap = union(swap, indexb);
%             swap = union(swap, local_orb(2));
            par = setdiff(par, indexb);
        else
            for innerCount = 1:N
                innerValue = B(innerCount, par(1));
                B(innerCount, par(1)) = B(innerCount, local_orb(2));
                B(innerCount, local_orb(2)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', par(1), local_orb(2));
            swap = union(swap, par(1));
%             swap = union(swap, local_orb(2));
            par = setdiff(par, par(1));
            % Assign the other particle indexa to local_orb(1)
            for innerCount = 1:N
                innerValue = B(innerCount, indexa);
                B(innerCount, indexa) = B(innerCount, local_orb(1));
                B(innerCount, local_orb(1)) = innerValue;
            end
            fprintf(datafile, 'Swap %d and %d \n', indexb, local_orb(1));
            swap = union(swap, indexa);
%             swap = union(swap, local_orb(1));
            par = setdiff(par, indexa);
        end
    elseif length(local_orb) > 2
        fprintf(datafile, 'Error! There are more than two orbitals inside the level 1 of the particle \n');
        break;
    elseif isempty(local_orb) == 1
        fprintf(datafile, 'Zero orbital inside the level 1 of the particle %d after cutoof \n  Then we use level 2 cutoff. \n \n', par(1));
        C = B;
        for inner = 1:length(swap)
            C(par(1), swap(inner)) = 0.0;
        end
        [maxValue, indexc] = max(abs(C(par(1), :)));
        % Assign indexc to par(1)
        for innerCount = 1:N
            innerValue = B(innerCount, indexc);
            B(innerCount, indexc) = B(innerCount, par(1));
            B(innerCount, par(1)) = innerValue;
        end
        fprintf(datafile, 'So, swap %d and %d \n', indexc, par(1));
        swap = union(swap, indexc);
        par = setdiff(par, par(1));
    end         
    
        [row, local_orb] = find(B(par(1), :) > ele_level_1);
    
fclose(datafile);
