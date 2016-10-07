global constK N sideln max_dist coe period

N = 175; dim = 1;
totalSteps = 50;
% n_size = 3:2:29;
% N = n_size;
constK = 0.5;
coe = 0.3;
[Pconfig_ideal, Oconfig, period, sigma, sideln, max_dist] = ...
    read_configuration_max_dis (1024, dim, N);

perturb = coe*sideln;
temp = Pconfig_ideal + perturb;
Pconfig_ideal_shifted = put_in_box(temp);

A_ideal_shifted = compute_matrix (Pconfig_ideal_shifted, Oconfig);
figure(1)
%     set(gca,'XTick',n_size)
%     ylim([-1 n_size(length(n_size))+1])
plot(eig(A_ideal_shifted), '*r');
% axis([-0.4 0.5 -0.4 0.4]);
hold on;
figure(totalSteps)
Pconfig = Pconfig_ideal_shifted;
E_nrm = zeros(length(totalSteps),1);
E_per = zeros(length(totalSteps),1);
cond_num = zeros(length(totalSteps),1);
% Get figure size
% pos = get(gcf, 'Position');
% width = pos(3); height = pos(4);

% Preallocate data (for storing frame data)
% imind = zeros(height, width, 1, length(totalSteps), 'uint8');

for i = 1:totalSteps
%     per_dis = unifrnd(-1,1, N, 1)*0.2*sideln;
%     perturb = normrnd(0, 0.1)*sideln;
    per_dis = normrnd(0, 1, N, 1)*0.5*sideln;
    temp = Pconfig + per_dis;
    Pconfig = put_in_box(temp);

    A = compute_matrix (Pconfig, Oconfig);
    E = A - A_ideal_shifted;
    E_nrm(i) = norm(E);
%     dlmwrite(sprintf('%dA.txt', i), A);
    cond_num(i) = cond(A);
%     if (i == 1) && (i == 2)
%     if i < 2
        plot(eig(A_ideal_shifted), '*r');
        hold on
        plot(eig(A), '*b');
        axis([-0.4 0.5 -0.4 0.4]);
        frame(i) = getframe(totalSteps);
        hold off
%     else
%         plot(eig(A_ideal_shifted), '*r');
%         hold on
%         plot(eig(A), '*b');
%         axis([-0.4 0.5 -0.4 0.4]);
%         frame = getframe(totalSteps);
%         hold off
%     end
%     im = frame2im(frame);
%     if i == 1
%         [imind(:,:,1,i),map] = rgb2ind(frame.cdata,256,'nodither');
% %         [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
%     else
%         imind(:,:,1,i) = rgb2ind(frame.cdata, map,'nodither');
% %         mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
%     end
    
end
movie2gif(frame, 'test.gif', 'LoopCount', 0, 'DelayTime', 0)
% imwrite(imind,map,'testgif.gif','DelayTime',0,'LoopCount',inf);
figure(totalSteps+1)
plot(log10(cond_num), '-*b');
figure(totalSteps+2)
plot(E_nrm, '-+b');