function [Pconfig, Oconfig, Box, sigma, constK, sideln, max_dist_sqr, max_Acut] = read_configuration (systemsize)

global dim

% 2012-04-18 add sidlength to configurate particles;
% 2012-11-07 add max_dist as a cutoff radius when computing A and u.

N = systemsize;
constK = 1.0;
% constK = 0.5;
%constK = 0.35;
max_Acut = 1e+6; % 1/(max_Acut) is the smallest nonzero value of A, like 1/1e+6 = 1e-6.
max_dist_sqr = log(max_Acut/3.5449024699805363)/constK;  % This is the cutoff radius when calculate A and u.
sideln = 2.030982595126518;

if dim == 3
    % Using BCC construction.
    nrK = floor((N/2)^(1.0/3.0)+1e-4); % Number of orbitals along one side of the box.
    Box = nrK*sideln*[1 1 1];
    sigma = Box(1)*0.02;
elseif dim == 2
    n = floor(sqrt(N)+1e-3);
    Box = (n*sideln)*[1 1];
    sigma = sideln*(0.5 - 0.1*dim);
elseif dim == 1
    Box = (N*sideln);
    sigma = sideln*(0.5 - 0.1*dim);
end

% Oconfig = buildBCC (N, sideln);
Oconfig = dlmread(sprintf('orbital%d',N));
% Oconfig = dlmread('K=1_15000_Orbital.txt');
Oconfig = Oconfig(:,1:3);
% Pconfig = dlmread(sprintf('particle%d',N));
Pconfig = dlmread('particle30K.txt');
Pconfig = Pconfig(:,1:3);
end


%----------------------------------------------------------
function Oconfig = buildBCC (N, sideln)

global dim

Oconfig = zeros(N, dim);
if dim == 1
    for j = 1:num
    Oconfig(j) = (j-1)*sideln;
    end
elseif dim == 2
    n = floor(sqrt(N)+1e-3);
    for j = 1:n^2
        r = zeros(1, 2);
        ip = floor((j-1)/n);
        r(1) = ip*sideln;
        r(2) = (j-1 - ip*n)*sideln;
        Oconfig(j, :) = r;
    end
elseif dim == 3
    n = floor(N^(1/3)+1e-3);
    for j = 1:n^3
        r = zeros(1, 3);
        ip = floor((j-1)/n^2);
        r(1) = ip*sideln;
        ip2 = floor((j-1 - ip*n^2)/n);
        r(2) = ip2*sideln;
        r(3) = (j-1 - ip*n^2 - ip2*n)*sideln;
        Oconfig(j, :) = r;
    end
end
end