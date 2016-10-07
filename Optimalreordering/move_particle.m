function [u, Pconfig1] = move_particle (A, movingP, move, Pconfig, Oconfig)

global constK max_dist dim max_Acut

Pconfig1 = Pconfig;
N = size(Pconfig,1);

if dim == 1
    Pconfig1(movingP) = put_in_box(Pconfig1(movingP) + move);
else
    Pconfig1(movingP,:) = put_in_box(Pconfig1(movingP,:) + move);
end

u = size(N,1);
for j=1:N
    if dim == 1
        d = convert_dist(Pconfig1(movingP) - Oconfig(j));
    else
        d = convert_dist(Pconfig1(movingP,:) - Oconfig(j,:));
    end
    dist = d*d';
    if (dist < max_dist)
        u(j, 1) = fix((1/3.5449024699805363) * exp(-constK*dist)*max_Acut)/max_Acut - A(movingP, j);
    else
        % 2012-12-18. change.
        u(j, 1) = 0.0 - A(movingP, j);
%         u(j, 1) = 0.0;
    end
%     clear d dist;
end

