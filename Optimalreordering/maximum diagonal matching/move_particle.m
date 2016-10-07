function [u, Pconfig1] = move_particle (A, movingP, move, Pconfig, Oconfig, dim)

global constK max_dist

Pconfig1 = Pconfig;
N = size(Pconfig,1);

if dim == 1
    Pconfig1(movingP) = put_in_box(Pconfig1(movingP) + move);
elseif dim == 2
    Pconfig1(movingP,:) = put_in_box(Pconfig1(movingP,:) + move);
elseif dim == 3
    Pconfig1(movingP,:) = put_in_box(Pconfig1(movingP,:) + move);
end
% Pconfig1(movingP) = put_in_box(Pconfig1(movingP) + move);

u = size(N,1);
for j=1:N
    d = convert_dist(Pconfig1(movingP,:) - Oconfig(j,:));
%     d = convert_dist(Pconfig1(movingP) - Oconfig(j));
    dist = d*d';
    if (dist < max_dist)
        u(j, 1) = (1/3.5449024699805363) * exp(-constK*dist) - A(movingP, j);
    else
        % 2012-12-18. change.
        u(j, 1) = 0.0 - A(movingP, j);
%         u(j, 1) = 0.0;
    end
%     clear d dist;
end

