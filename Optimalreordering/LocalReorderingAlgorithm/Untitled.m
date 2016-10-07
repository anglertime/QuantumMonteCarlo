
global constK max_dist_sqr dim max_Acut
N = 1024;
M = A5;
partC = [ 12.034449 14.594160 10.091691 ];
u = size(N,1);
for j=1:N
    if dim == 1
        d = convert_dist(partC - OrbitalConfig(j));
    else
        d = convert_dist(partC - OrbitalConfig(j,:));
    end
    dist = d*d';
    if (dist < max_dist_sqr)
        u(j, 1) = fix((1/3.5449024699805363) * exp(-constK*dist)*max_Acut)/max_Acut - M(891, j);
    else
        % 2012-12-18. change.
        u(j, 1) = 0.0 - M(891, j);
%         u(j, 1) = 0.0;
    end
%     clear d dist;
end

