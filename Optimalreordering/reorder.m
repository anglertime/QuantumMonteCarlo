global constK N_length sideln max_dist coe period

N_length = 40; dim = 1;
totalSteps = 1;
% n_size = 3:2:29;
% N = n_size;
constK = 0.5;
coe = 0.3;
[Pconfig_ideal, Oconfig, period, sigma, sideln, max_dist] = ...
    read_configuration_max_dis (1024, dim, N_length);

perturb = coe*sideln;
temp = Pconfig_ideal + perturb;
Pconfig_ideal = put_in_box(temp);

Pconfig_ideal_moved = Pconfig_ideal;
Pconfig_ideal_moved(25) = Pconfig_ideal_moved(35) - 25.1*sideln;
Pconfig_ideal_moved = put_in_box(Pconfig_ideal_moved);
% Pconfig_ideal_moved(26) = Pconfig_ideal_moved(26) - 14.7*sideln;

A_ideal = compute_matrix (Pconfig_ideal, Oconfig);
A_ideal_moved = compute_matrix (Pconfig_ideal_moved, Oconfig);
figure(1)
%     set(gca,'XTick',n_size)
%     ylim([-1 n_size(length(n_size))+1])
plot(real(eig(A_ideal)),imag(eig(A_ideal)), '*r')
axis([-0.3 0.4 -0.2 0.2]);
hold on
plot(real(eig(A_ideal_moved)),imag(eig(A_ideal_moved)), 'ob')
% axis([-0.4 0.5 -0.4 0.4]);


% S = svds(A_ideal_moved,40);
% S(40) 
% S(1)
cond(A_ideal_moved)
[A_do] = make_diagonal_dominant(A_ideal_moved);
% [A_re] = reorder_try_min_maxdis(A_ideal_moved);
[A_shift_re] = shift_reorder( A_ideal_moved );

figure(2); 

plot(real(eig(A_ideal_moved)),imag(eig(A_ideal_moved)), 'ob');
axis([-0.3 0.4 -0.2 0.2]);
hold on; plot(real(eig(A_do)),imag(eig(A_do)), '+r'); 
hold on; plot(real(eig(A_shift_re)),imag(eig(A_shift_re)), '*k')
% figure(totalSteps)
% Pconfig = Pconfig_ideal;
% E_nrm = zeros(length(totalSteps),1);
% E_per = zeros(length(totalSteps),1);
% cond_num = zeros(length(totalSteps),1);
% % Get figure size
% % pos = get(gcf, 'Position');
% % width = pos(3); height = pos(4);
% 
% % Preallocate data (for storing frame data)
% % imind = zeros(height, width, 1, length(totalSteps), 'uint8');
% 
% for i = 1:N_length
% %     temp = Pconfig_ideal(i) + perturb;
% %     Pconfig_ideal_shifted = put_in_box(temp);
% %     per_dis = unifrnd(-1,1, N_length, 1)*0.2*sideln;
% %     perturb = normrnd(0, 0.1)*sideln;
% %     per_dis = normrnd(0, 0.1, N, 1)*coe*sideln;
%     temp = Pconfig(i) + perturb;
%     Pconfig(i) = put_in_box(temp);
% 
%     A = compute_matrix (Pconfig, Oconfig);
%     E = A - A_ideal_shifted;
%     E_nrm(i) = norm(E);
% %     dlmwrite(sprintf('%dA.txt', i), A);
%     cond_num(i) = cond(A);
% %     if (i == 1) && (i == 2)
% %     if i < 2
%         plot(eig(A_ideal_shifted), '*r');
%         hold on
%         plot(eig(A), '*b');
%         axis([-0.4 0.5 -0.4 0.4]);
%         frame(i) = getframe(totalSteps);
%         s(i) = svds(A,1,0);
%         nrm_1(i) = norm(A, 1);
%         nrm_inf(i) = norm(A, Inf);
%         hold off
% %     else
% %         plot(eig(A_ideal_shifted), '*r');
% %         hold on
% %         plot(eig(A), '*b');
% %         axis([-0.4 0.5 -0.4 0.4]);
% %         frame = getframe(totalSteps);
% %         hold off
% %     end
% %     im = frame2im(frame);
% %     if i == 1
% %         [imind(:,:,1,i),map] = rgb2ind(frame.cdata,256,'nodither');
% % %         [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
% %     else
% %         imind(:,:,1,i) = rgb2ind(frame.cdata, map,'nodither');
% % %         mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
% %     end
%     
% end
% s = s(1:N_length);
% nrm_1 = nrm_1(1:N_length);
% nrm_inf = nrm_inf(1:N_length);
% 
% figure(N_length+1)
% plot(log10(cond_num), '-+b')
% hold on;
% title({sprintf('cond num/E norm. k=0.5 %dD n=%d', dim, N_length)});
% plot(log10(E_nrm), '-*r')
% hold on
% xlabel({'step'});
% ylabel({'log10 magnitude'});
% 
% figure(N_length+2)
% plot((nrm_1), '-+b')
% title({sprintf('1-norm/Inf-norm. k=0.5 %dD n=%d', dim, N_length)});
% hold on;
% plot((nrm_inf), '-*r')
% hold on
% xlabel({'step'});
% ylabel({'value'});
% 
% movie2gif(frame, 'test.gif', 'LoopCount', 0, 'DelayTime', 0)
% % imwrite(imind,map,'testgif.gif','DelayTime',0,'LoopCount',inf);
% figure(length(totalSteps)+1)
% plot(cond_num, '-*b');
% figure(length(totalSteps)+2)
% plot(E_nrm, '-+b');