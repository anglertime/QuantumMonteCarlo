function [ B ] = shift_reorder( B )

global sideln sideLength constK N_length max_dist period
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

datafile = fopen(sprintf('shift reorder K%d.txt',constK),'wt+');

level_1 = 0.50001*sideln;
% level_1 = 0.5001*sideLength;
ele_level_1 = (1/3.5449024699805363) * exp(-constK*level_1^2);
orb_nbr = zeros(N_length, 1);
k = 1;
for i = 1:N_length
    [row, local_orb] = find(B(:, i) > ele_level_1);
    orb_nbr(i) = length(row);
    if orb_nbr(i) > 2
        fprintf(datafile, 'Error! There are more than two particles around orbital %d\n', i);
        break;
    elseif orb_nbr(i) == 2
        orb = i;
        par1 = row(1);
        par2 = row(2);
        fprintf(datafile, 'Orb is %d. Par1 is %d and par2 is %d. \n', orb, par1, par2);
        k = k + 1;
        if k > 3
            fprintf(datafile, 'k > 3. There are two orbitals having two particles \n');
            break;
        end
    elseif orb_nbr(i) == 0
        par = i;
        fprintf(datafile, 'zero orb %d \n', i);
    end
end

fprintf(datafile, 'orb is %d. par is %d \n', orb, par);
% if mod(abs(orb - par), N_length) == 1
%     
% end
if abs(orb - par) < floor(N_length/2 + 1.001)
    if par > orb
        temp = B(par, :);
        for innerCount = par:-1:orb+2
            B(innerCount, :) = B(innerCount - 1, :);
%             fprintf(datafile, 'Swap. %d \n', innerCount);
        end
        B(orb+1, :) = temp;
    elseif orb > par
        temp = B(par, :);
        for innerCount = par:orb-2
            B(innerCount, :) = B(innerCount + 1, :);
        end
        B(orb-1, :) = temp;
    end
elseif abs(orb - par) > floor(N_length/2 + 1.001)
    if par > orb
        temp1 = B(par, :);
        temp2 = B(1, :);
        for innerCount = 1:orb-2
            B(innerCount, :) = B(innerCount + 1, :);
        end
        B(orb-1, :) = temp1;
        for innerCount = par:N_length-1
            B(innerCount, :) = B(innerCount + 1, :);
        end
        B(N_length,:) = temp2;
    elseif orb > par
        temp1 = B(par, :);
%         temp2 = B(1, :);
        for innerCount = par:-1:2
            B(innerCount, :) = B(innerCount - 1, :);
        end
        B(1, :) = B(N_length, :);
        for innerCount = N_length:-1orb+2
            B(innerCount, :) = B(innerCount - 1, :);
        end
        B(orb+1,:) = temp1;
    end
end
    


fclose(datafile);
end

