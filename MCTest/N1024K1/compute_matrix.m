function  A = compute_matrix (Pconfig, Oconfig)

global constK max_dist_sqr dim max_Acut

if dim == 1
    N = length(Pconfig);
else
    N = size(Pconfig, 1);
end

k = 1;
for i=1:N
    for j=1:N
        d = convert_dist(Pconfig(i, :) - Oconfig(j, :));
        dist = d*d';
        if (dist < max_dist_sqr)
            c(k) = i;
            r(k) = j;
            % Round the value into 1/max_Acut
            v(k) = fix((1/3.5449024699805363) * exp(-constK*dist)*max_Acut)/max_Acut;
            k = k+1;
        end
    end
end
A = sparse(c, r, v);
A = full(A);
end