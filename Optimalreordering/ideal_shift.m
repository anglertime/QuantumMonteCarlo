global constK N sideln max_dist period

N = 225; dim = 2;

per_dis = 0:0.1:2.6;
% n_size = 3:2:29;
% N = n_size;
constK = 0.5;
[Pconfig_ideal, Oconfig, period, sigma, sideln, max_dist] = ...
    read_configuration_max_dis (1024, dim, N);

A_ideal = compute_matrix (Pconfig_ideal, Oconfig);
figure(1)
%     set(gca,'XTick',n_size)
%     ylim([-1 n_size(length(n_size))+1])
plot(eig(A_ideal), '*r');
% axis([-0.4 0.5 -0.4 0.4]);
hold on;
figure(length(per_dis))
cond_num = zeros(length(per_dis),1);
% Get figure size
% pos = get(gcf, 'Position');
% width = pos(3); height = pos(4);

% Preallocate data (for storing frame data)
% imind = zeros(height, width, 1, length(totalSteps), 'uint8');

for i = 1:length(per_dis)
    
    perturb = per_dis(i)*sideln;
%     Pconfig = Pconfig_ideal;
    Pconfig = Pconfig_ideal + sqrt(2)*0.5*perturb;
%     Pconfig = Pconfig_ideal + sqrt(3)*perturb/3;
%     Pconfig(:, 1) = Pconfig_ideal(:, 1) + perturb;
%     Pconfig(:, 2) = Pconfig_ideal(:, 2) + sqrt(2)*0.5*perturb;
%     Pconfig(:, 1) = Pconfig_ideal(:, 1) + sqrt(2)*0.5*perturb;
    
    Pconfig = put_in_box(Pconfig);

    A = compute_matrix (Pconfig, Oconfig);
%     dlmwrite(sprintf('%dA.txt', i), A);
    cond_num(i) = cond(A);
%     if (i == 1) && (i == 2)
%     if i < 2
        plot(eig(A_ideal), '--*r');
        hold on
        plot(eig(A), '*b');
        axis([-0.4 0.5 -0.4 0.4]);
        frame(i) = getframe(length(per_dis));
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
figure(length(per_dis)+1)
plot(log10(cond_num), '-*b');