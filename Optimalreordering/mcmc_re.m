global constK datafile detRatio_dense sideLength max_dist period

size = 125; dim = 1;
constK = 0.5;
[ParticleConfig, OrbitalConfig, period, sigma, sideLength, max_dist] = ...
    read_configuration_max_dis (1024, dim, size);
per_dis = normrnd(0, 1, size, dim)*0.2*sideLength;
per_dis = per_dis + 0.2*sideLength;
temp = ParticleConfig + per_dis;
ParticleConfig = put_in_box(temp);

% [ParticleConfig, OrbitalConfig, Box, sigma, constK, sideLength, max_dist] = read_configuration;

totalSteps = 300;
datafile = fopen(sprintf('MCMC Test Dim%d K%d.txt', dim, constK),'wt+');

% N = size(ParticleConfig,1);
N = length(ParticleConfig);
A = compute_matrix (ParticleConfig, OrbitalConfig);
A_ideal = A;
B = A;
% A = dlmread('12000Matrix.txt');

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

pass = 0;
count_each_part(N)=0;
acceptCount = 0;
    
[movingParticle, move] = random_move (N, sigma, dim);
randNum = rand;

% nr_v(totalSteps) = 0;
% nr_in(totalSteps) = 0;
% test_result(totalSteps) = 0;
% error(totalSteps) = 0;
% error1(totalSteps) = 0;
figure(totalSteps)
kk = 1;
for mcStep = 1:totalSteps
    


   [u, ParticleConfig1] = move_particle (A, movingParticle, move, ParticleConfig, OrbitalConfig, dim);

   fprintf(datafile, 'At step %d \n \tmoving part %d from %s to %s \n', mcStep,...
       movingParticle, pos_to_string(ParticleConfig(movingParticle, :), dim), pos_to_string(ParticleConfig1(movingParticle, :), dim));

   %--------------  Dense Method  --------------------
   detRatio_dense = dense_detRatio(A,u,movingParticle);

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

passb = 0;
count_each_partb(N)=0;
acceptCountb = 0;
figure(2*totalSteps)
kkk = 1;
for mcStep = 1:totalSteps
    


   [u, ParticleConfig1] = move_particle (B, movingParticle, move, ParticleConfig, OrbitalConfig, dim);

   fprintf(datafile, 'At step %d \n \tmoving part %d from %s to %s \n', mcStep,...
       movingParticle, pos_to_string(ParticleConfig(movingParticle, :), dim), pos_to_string(ParticleConfig1(movingParticle, :), dim));

   %--------------  Dense Method  --------------------
   detRatio_dense = dense_detRatio(B,u,movingParticle);

   detRatio = detRatio_dense;
   if (detRatio*detRatio > randNum)
       fprintf(datafile, '-----accept   %f   %f \n\n', detRatio*detRatio, randNum);
       B(movingParticle,:) = B(movingParticle,:) + u';
       ParticleConfig = ParticleConfig1;
       acceptCountb = acceptCountb+1;
        plot(eig(B_ideal), '*r');
        hold on
        plot(eig(B), '*b');
        axis([-0.4 0.5 -0.4 0.4]);
        frameb(kkk) = getframe(2*totalSteps);
        cond_nrb(kkk) = cond(B);
        sb(kkk) = svds(B,1,0);
        nrm_1b(kkk) = norm(B, 1);
        nrm_infb(kkk) = norm(B, Inf);
        hold off
        kkk = kkk + 1;
        [ B, ~ ] = max_min_diagonal_final( B );
   else
       fprintf(datafile, '-----reject    %f   %f \n\n', detRatio*detRatio, randNum);
   end

    % Find if a pass is completed
    for count = 1:N
        count_each_partb(count) = length(find(movingParticleArr == count));
    end

    % Check the observables every pass exactly
    if (rem(mcStep,N)== 0 && min(count_each_partb) > passb)
        passb = passb + 1;
        %CalcPairCorrelation(ParticleConfig, pairHistogram);
    end
end
cond_nrb = cond_nrb(1:acceptCountb);
sb = sb(1:acceptCountb);
nrm_1b = nrm_1b(1:acceptCountb);
nrm_infb = nrm_infb(1:acceptCountb);

movie2gif(frameb, sprintf('k=0.5 %dD n=%d shifted testb.gif', dim, size), 'LoopCount', 0, 'DelayTime', 0)
% movie2gif(frame, 'test.gif', 'LoopCount', 0, 'DelayTime', 0)
figure(2*totalSteps+1)
plot(log10(cond_nrb), '-+b')
hold on;
title({sprintf('cond num/smallest sing val. k=0.5 %dD n=%d', dim, size)});
plot(log10(sb), '-*r')
hold on
plot(log10(cond_nrb.*sb), '--+m')
xlabel({'Accepted MC step'});
ylabel({'log10 magnitude'});


figure(2*totalSteps+2)
plot((nrm_1b), '-+b')
title({sprintf('1-norm/Inf-norm. k=0.5 %dD n=%d', dim, size)});
hold on;
plot((nrm_infb), '-*r')
hold on
xlabel({'Accepted MC step'});
ylabel({'value'});

fprintf(datafile, '\n %d out of %d movements are accepted. Accept ratio is %6.4f%%  \n', ...
    acceptCountb, totalSteps, 100*acceptCountb/totalSteps);
fprintf(datafile, '\n Finally, after %d steps, pass = %d. \n', totalSteps, passb);
           
fclose(datafile);