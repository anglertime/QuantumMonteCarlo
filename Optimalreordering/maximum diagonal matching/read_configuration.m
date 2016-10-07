function [Pconfig, Oconfig, Box, sigma, constK, sideln, max_dist] = read_configuration ()

% 2012-04-18 add sidlength to configurate particles;
% 2012-11-07 add max_dist as a cutoff radius when computing A and u.

N = 2*6^3;
% N = 1024;
%N = 1458;
%N = 2000;
nrK = floor((N/2)^(1.0/3.0)+1e-4); % Number of orbitals along one side of the box.

constK = 1.0;
% constK = 0.5;
%constK = 0.35;
max_dist = 5e+5; % 1/(max_dist) is the smallest nonzero value of A, like 1/5e+5 = 2e-6.
max_dist = log(max_dist/3.5449024699805363)/constK;  % This is the cutoff radius when calculate A and u.

Box = (N*4.0/3.0*pi)^(1.0/3.0) * [1 1 1];
sideln = Box(1)/nrK;
sigma = Box(1)*0.02;

Oconfig = buildBCC (N, Box);
% Oconfig = dlmread(sprintf('orbital%d',N));
% Oconfig = dlmread('K=1_15000_Orbital.txt');
Oconfig = Oconfig(:,1:3);
% Pconfig = dlmread(sprintf('particle%d',N));
% Pconfig = dlmread('30300particle_global.txt');

perturb = unidrnd(100, size(Oconfig))*0.1*sideln*0.01;
Pconfig = Oconfig + perturb;  % Start from the same configuration.
%for i = 1:N
%    move = randn(1,3)*sigma;
%    Pconfig(i, :) = Pconfig(i, :) + move;
%end
Pconfig = Pconfig(:,1:3);


%----------------------------------------------------------
function config = buildBCC (N, Box)

config = zeros(N, 3);
%NperDim = ceil ( (0.5*N)^(1.0/3.0) - 1e-6 );
% 2012-05-14 I change it
NperDim = floor ( (0.5*N)^(1.0/3.0) + 1e-4 );
delta = Box(1) / NperDim;

x = zeros(1,3);
for i=1:N
    ip = floor((i-1)/2);
    x(1) = floor(ip/(NperDim*NperDim));
    x(2) = floor((ip-x(1)*NperDim*NperDim)/NperDim);
    x(3) = floor(ip - x(1)*NperDim*NperDim - x(2)*NperDim);
    % 2012-05-14 delete '- 0.5*Box'
    r = x*delta;
    if (mod(i,2)==0)
        r = r + 0.5*delta;
    end
    config(i,:) = r;
end