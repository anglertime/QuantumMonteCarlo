function [Pconfig_ideal, Oconfig, period, sigma, sideln, max_dist] = read_configuration_max_dis (sys_size, dim, size)

global constK

nrK = floor((sys_size/2)^(1.0/3.0)+1e-4); % Number of orbitals along one side of the box.

max_dist = 5e+5; % 1/(max_dist) is the smallest nonzero value of A, like 1/5e+5 = 2e-6.
max_dist = log(max_dist/3.5449024699805363)/constK;  % This is the cutoff radius when calculate A and u.
Box = (sys_size*4.0/3.0*pi)^(1.0/3.0) * [1 1 1];
sideln = Box(1)/nrK;
% sigma = Box(1)*0.02;

num = size;
if dim == 1
    Oconfig = 0:sideln:(num-1)*sideln;
    Oconfig = Oconfig';
    Pconfig_ideal = Oconfig;  % Start from the same configuration.
    period = num*sideln;
elseif dim == 2
    num = sqrt(num);
    % 2D
    Oconfig = zeros(num^2, 2);
    for j = 1:num^2
        r = zeros(1, 2);
        ip = floor((j-1)/num);
        r(1) = ip*sideln;
        r(2) = (j-1 - ip*num)*sideln;
        Oconfig(j, :) = r;
    end
    Pconfig_ideal = Oconfig;  % Start from the same configuration.
    period = num*sideln;
elseif dim == 3
    num = floor(num^(1/3)+0.01);
    % 3D
    Oconfig = zeros(num^3, 3);
    for j = 1:num^3
        r = zeros(1, 3);
        ip = floor((j-1)/num^2);
        r(1) = ip*sideln;
        ip2 = floor((j-1 - ip*num^2)/num);
        r(2) = ip2*sideln;
        r(3) = (j-1 - ip*num^2 - ip2*num)*sideln;
        Oconfig(j, :) = r;
    end
    Pconfig_ideal = Oconfig;  % Start from the same configuration.
    period = num*sideln;
else
end
sigma = 0.6*sideln;