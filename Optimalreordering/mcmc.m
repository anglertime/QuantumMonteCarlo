global constK datafile detRatio_dense sideLength max_dist period

size = 115; dim = 1;
constK = 0.5;
[ParticleConfig, OrbitalConfig, period, sigma, sideLength, max_dist] = ...
    read_configuration_max_dis (1024, dim, size);
per_dis = normrnd(0, 1, size, dim)*0.2*sideLength;
per_dis = per_dis + 0.2*sideLength;
temp = ParticleConfig + per_dis;
ParticleConfig = put_in_box(temp);

% [ParticleConfig, OrbitalConfig, Box, sigma, constK, sideLength, max_dist] = read_configuration;

totalSteps = 500;
datafile = fopen(sprintf('MCMC Test Dim%d K%d.txt', dim, constK),'wt+');

% N = size(ParticleConfig,1);
N = length(ParticleConfig);
A = compute_matrix (ParticleConfig, OrbitalConfig);
A_ideal = A;
B = A;
% A = dlmread('12000Matrix.txt');

disp(sprintf('Max element in first row = %f \n', max(A(1,:))));
disp(sprintf('Max element in first column = %f \n', max(A(:,1))));
for tempCount = 1:10
    disp(sprintf('%e ', A(1, tempCount)));
end
disp(sprintf('\n'));
for tempCount = 1:10
    disp(sprintf('%e ', A(tempCount, 1)));
end
maxN = 3; 
pairHistogram(1:100) = 0; 
totalHistogram(1:maxN, 1:maxN, 1:maxN) = 0; 
totalEnergy = 0;

%This is for debugging that why energy and sf for the 300 system are
%not correct.

%Box
% boxData = fopen('Box.txt','wt+');
% fprintf(boxData, '%f %f %f', period, period, period); 
% fclose(boxData);        
%Initial coordinates
% initialStep = fopen('initialParticle.txt','wt+');
% save_config(initialStep, N, ParticleConfig);
% fclose(initialStep);    

movingParticleArr(totalSteps)=0;
pass = 0;
count_each_part(N)=0;
acceptCount = 0;

% nr_v(totalSteps) = 0;
% nr_in(totalSteps) = 0;
% test_result(totalSteps) = 0;
% error(totalSteps) = 0;
% error1(totalSteps) = 0;
figure(totalSteps)
kk = 1;
for mcStep = 1:totalSteps
    
    
   [movingParticle, move] = random_move (N, sigma, dim);
   randNum = rand;
   movingParticleArr(mcStep) = movingParticle;

   [u, ParticleConfig1] = move_particle (A, movingParticle, move, ParticleConfig, OrbitalConfig, dim);

   fprintf(datafile, 'At step %d \n \tmoving part %d from %s to %s \n', mcStep,...
       movingParticle, pos_to_string(ParticleConfig(movingParticle, :), dim), pos_to_string(ParticleConfig1(movingParticle, :), dim));

   %--------------  Dense Method  --------------------
   detRatio_dense = dense_detRatio(A,u,movingParticle);
   
   %--------------  Iterative Method  -------------------
%    [A11, A12, A21, A22, detRatio_iter, iter, err, err1, flag, r_nrm, innertol, inner_nr, in_iter] = ine_stgmres(1e-3, A, u, movingParticle, 1e-4);
%    [detRatio_iter, iter, err, residual, flag, r_nrm, v1pn, vpn] = gmres_sol(A, u, mcStep, movingParticle);
%     
%     AA = A;
%     nr_v(mcStep) = iter;
%     nr_in(mcStep) = in_iter;
%     error(mcStep) = err;
%     error1(mcStep) = err1;
%     f1 = detRatio_dense*detRatio_dense;
%     f2 = detRatio_iter*detRatio_iter;
%     f = abs(min(f1, 1) - min(f2, 1));
%     test_result(mcStep) = f;
%     if (err > 5e-3)
%         fprintf(datafile1, 'At step %d, error is %8.4e and f is %8.4e  \n', mcStep, err, f);
%     end
% 
%    fprintf(datafile, '\tdense: %f   iterative:%f \n', detRatio_dense, detRatio_iter);

   detRatio = detRatio_dense;
   if (detRatio*detRatio > randNum)
       fprintf(datafile, '-----accept   %f   %f \n\n', detRatio*detRatio, randNum);
       A(movingParticle,:) = A(movingParticle,:) + u';
       ParticleConfig = ParticleConfig1;
       acceptCount = acceptCount+1;
        plot(eig(A_ideal), '*r');
        hold on
        plot(eig(A), '*b');
        axis([-0.4 0.5 -0.4 0.4]);
        frame(kk) = getframe(totalSteps);
        cond_nr(kk) = cond(A);
        s(kk) = svds(A,1,0);
        nrm_1(kk) = norm(A, 1);
        nrm_inf(kk) = norm(A, Inf);
        hold off
        kk = kk + 1;
   else
       fprintf(datafile, '-----reject    %f   %f \n\n', detRatio*detRatio, randNum);
   end

    % Find if a pass is completed
    for count = 1:N
        count_each_part(count) = length(find(movingParticleArr == count));
    end

    % Check the observables every pass exactly
    if (rem(mcStep,N)== 0 && min(count_each_part) > pass)
        pass = pass + 1;
        %CalcPairCorrelation(ParticleConfig, pairHistogram);
    end
    
    %Save configuration to check equilibrium
     if (mod(mcStep, 2000) == 0)
         Steps = fopen(sprintf('%dparticle.txt', mcStep), 'wt+');
         save_config(Steps, N, ParticleConfig);
         fclose(Steps);   
%          dlmwrite(sprintf('%derror.txt', mcStep), error);
%          dlmwrite(sprintf('%dnr_v.txt', mcStep), nr_v(1:mcStep));
%          dlmwrite(sprintf('%derror1.txt', mcStep), error1);
%          dlmwrite(sprintf('%dtest.txt', mcStep), test_result(1:mcStep));
         dlmwrite(sprintf('%dA.txt', mcStep), A);
         fprintf(datafile, sprintf('\n At step %d, accept move is %d, \n Accept ratio is %6.4f%%. \n', ...
             mcStep, acceptCount, 100*acceptCount/totalSteps));
     end
end
cond_nr = cond_nr(1:acceptCount);
s = s(1:acceptCount);
nrm_1 = nrm_1(1:acceptCount);
nrm_inf = nrm_inf(1:acceptCount);

movie2gif(frame, sprintf('k=0.5 %dD n=%d shifted test.gif', dim, size), 'LoopCount', 0, 'DelayTime', 0)
% movie2gif(frame, 'test.gif', 'LoopCount', 0, 'DelayTime', 0)
figure(totalSteps+1)
plot(log10(cond_nr), '-+b')
hold on;
title({sprintf('cond num/smallest sing val. k=0.5 %dD n=%d', dim, size)});
plot(log10(s), '-*r')
hold on
plot(log10(cond_nr.*s), '--+m')
xlabel({'Accepted MC step'});
ylabel({'log10 magnitude'});


figure(totalSteps+2)
plot((nrm_1), '-+b')
title({sprintf('1-norm/Inf-norm. k=0.5 %dD n=%d', dim, size)});
hold on;
plot((nrm_inf), '-*r')
hold on
xlabel({'Accepted MC step'});
ylabel({'value'});

fprintf(datafile, '\n %d out of %d movements are accepted. Accept ratio is %6.4f%%  \n', ...
    acceptCount, totalSteps, 100*acceptCount/totalSteps);
fprintf(datafile, '\n Finally, after %d steps, pass = %d. \n', totalSteps, pass);

% fifthSteps = fopen('FinalParticle.txt','wt+');
% save_config(fifthSteps, N, ParticleConfig);
% fclose(fifthSteps);            
fclose(datafile);