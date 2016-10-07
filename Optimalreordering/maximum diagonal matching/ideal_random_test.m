global constK N sideln max_dist coe

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
n_size = 90;
% per_dis = [0.1, 0.2, 0.4, 0.5, 0.6];
coe = 0.3;

% per_dis = normrnd(0, 0.1, 100, 1)*coe*sideln;
per_dis = unifrnd(-1,1, 100, 1)*0.15*coe*sideln;

N = n_size(1);

    Oconfig = 0:sideln:(N-1)*sideln;
    Oconfig = Oconfig';
    Pconfig_ideal = Oconfig;  % Start from the same configuration.
    perturb = coe*sideln;
    temp = Pconfig_ideal + perturb;
%     Pconfig_ideal_shifted = zeros(size(Pconfig_ideal));
    
    %1D
    Pconfig_ideal_shifted = put_in_box(temp);
    
%     for i = 1:N
%         Pconfig_ideal_shifted(i) = put_in_box(temp(i));
%     end
    
    A_ideal_shifted = compute_matrix (Pconfig_ideal_shifted, Oconfig);
    figure(length(per_dis)+1)
%     set(gca,'XTick',n_size)
%     ylim([-1 n_size(length(n_size))+1])
    plot(eig(A_ideal_shifted), '*r');
    hold on;
    figure(1)
%     filename = 'testgif.gif';

% Pconfig = Pconfig_ideal_shifted;
for i = 1:length(per_dis)
    

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
%     Oconfig = zeros(N^3, 3);
%     for j = 1:N^3
%         r = zeros(1, 3);
%         ip = floor((j-1)/N^2);
%         r(1) = ip*sideln;
%         ip2 = floor((j-1 - ip*N^2)/N);
%         r(2) = ip2*sideln;
%         r(3) = (j-1 - ip*N^2 - ip2*N)*sideln;
%         Oconfig(j, :) = r;
%     end
%     Pconfig = Oconfig;  % Start from the same configuration.
    
%     perturb = floor(unidrnd(99, length(Pconfig), 1) - 50)*0.02*per_dis(2)*sideln;
%     perturb = normrnd(0, 0.1)*sideln;
    temp = Pconfig_ideal_shifted + per_dis(i);
    Pconfig = put_in_box(temp);
    %1D
%     for j = 1:N
%         Pconfig(j,:) = put_in_box(temp(j,:));
%     end

    A = compute_matrix (Pconfig, Oconfig);
    E = A - A_ideal_shifted;
    E_nrm(i) = norm(E);
%     dlmwrite(sprintf('%dA.txt', i), A);
    cond_num(i) = cond(A);
%     if (i == 1) && (i == 2)
    if i < 10
        plot(eig(A_ideal_shifted), '*r');
        hold on
        plot(eig(A), '*b');
        axis([0.1 0.4 -0.12 0.12]);
        frame = getframe(1);
        hold off
    else
        plot(eig(A_ideal_shifted), '*r');
        hold on
        plot(eig(A), '*b');
        axis([0.1 0.4 -0.12 0.12]);
        frame = getframe(1);
        hold off
    end
%     im = frame2im(frame);
    if i == 1
        [imind(:,:,1,i),map] = rgb2ind(frame.cdata,256,'nodither');
%         [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        imind(:,:,1,i) = rgb2ind(frame.cdata, map,'nodither');
%         mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
    end
%     [imind,cm] = rgb2ind(frame.cdata,256,'nodither');
%     imwrite(imind,cm,'testgif.gif','DelayTime',1,'LoopCount',inf);
    
end
imwrite(imind,map,'testgif.gif','DelayTime',0,'LoopCount',inf);
figure(length(per_dis)+2)
% plot(per_dis, cond_num, '-ob');
plot(E_nrm, '-ob');