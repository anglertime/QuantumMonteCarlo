global constK N sideln max_dist

side_N = 1024;
nrK = floor((side_N/2)^(1.0/3.0)+1e-4); % Number of orbitals along one side of the box.

con = [0.35, 0.5, 0.75, 1];
constK = 0.5;
% constK = 1;
max_dist = 5e+5; % 1/(max_dist) is the smallest nonzero value of A, like 1/5e+5 = 2e-6.
max_dist = log(max_dist/3.5449024699805363)/constK;  % This is the cutoff radius when calculate A and u.
Box = (side_N*4.0/3.0*pi)^(1.0/3.0) * [1 1 1];
sideln = Box(1)/nrK;
sigma = Box(1)*0.02;
 
% n_size = [5, 10, 20, 50, 70, 100, 150, 200, 300];
% n_size = 3:2:29;
n_size = 2:5;
per_dis = [0.1, 0.2, 0.4, 0.5, 0.6];
perturb = floor(unidrnd(99, length(Pconfig), 1) - 50)*0.02*per_dis(2)*sideln;
% perturb = 0.2*sideln;
Pconfig = Pconfig + perturb;

for i = 1:length(n_size)
    N = n_size(i);
    
%     Oconfig = 0:sideln:(N-1)*sideln;
%     Oconfig = Oconfig';
%     Pconfig = Oconfig;  % Start from the same configuration.

% 2D
%     Oconfig = zeros(N^2, 2);
%     for j = 1:N^2
%         r = zeros(1, 2);
%         ip = floor((j-1)/N);
%         r(1) = ip*sideln;
%         r(2) = (j-1 - ip*N)*sideln;
%         Oconfig(j, :) = r;
%     end
%     Pconfig = Oconfig;  % Start from the same configuration.

% 3D
    Oconfig = zeros(N^3, 3);
    for j = 1:N^3
        r = zeros(1, 3);
        ip = floor((j-1)/N^2);
        r(1) = ip*sideln;
        ip2 = floor((j-1 - ip*N^2)/N);
        r(2) = ip2*sideln;
        r(3) = (j-1 - ip*N^2 - ip2*N)*sideln;
        Oconfig(j, :) = r;
    end
    Pconfig = Oconfig;  % Start from the same configuration.
    
    A = compute_matrix (Pconfig, Oconfig);
%     dlmwrite(sprintf('%dA.txt', i), A);
    cond_num(i) = cond(A);
%     figure(i)
%     plot(eig(A), '*r');
%     hold on;
    figure(length(n_size)+1)
%     set(gca,'XTick',n_size)
%     ylim([-1 n_size(length(n_size))+1])
    plot(n_size(i), eig(A), '-*r');
    hold on;
    
end
figure(length(n_size)+2)
plot(n_size, cond_num, '-ob');