% 
% load('C:\Downloads\MATLAB\DataConfiguration\Equi3D2000BCC20K.mat');
% load('C:\Downloads\MATLAB\DataConfiguration\Equi3D686BCC18sweeps.mat');
load('C:\Downloads\MATLAB\DataConfiguration\Equi3DBCC12swp2662.mat');

global constK Box datafile detRatio_dense sideLength max_dist_sqr max_Acut
global indicator dim reinitialize datafile1 mat_L mat_U
global array_m array_u array_detRatio globalCount

%Initialization
globalCount = 0; love = 0;
dim = 3; N = 2*11^3; indicator = ones(1, N); totalSteps = 2*N;
array_m = [];  array_u = [];  array_detRatio = [];
sigma = 0.135*sideLength;
% ParticleConfig = Pconfig; OrbitalConfig = Oconfig;

% [ParticleConfig, OrbitalConfig, Box, sigma, constK, sideLength, max_dist_sqr, max_Acut] = read_configuration (N);
% A = compute_matrix (ParticleConfig, OrbitalConfig);
% Ainitial = A;
% A = Ainitial;
% load('30000A', 'A');
%% Initialization
% [ global_cutoff, ~ ] = max_min_diagonal ( A, 2*N+2, 5e-3, 1e-3);
% global_cutoff = 0.015;
localcutoff = global_cutoff;
% [ A, ParticleConfig, matching ] = max_weighted_match( A, cutoff, 3, ParticleConfig );
% dmfdfin = min(diag(A));
reinitialize = 0;
dsum = sum(diag(A));
%dsumcut = 1; % 1.2 > 1.
% if N < 500
%     cycle = N;
% else
%     cycle = floor(0.8*N);
% end

datafile = fopen(sprintf('Reordering vmc-OutputK%d.txt',constK),'wt+');
datafile1 = fopen(sprintf('Reordering vmc Large-Error-monitor-OutputK%d.txt',constK),'wt+');
disp(sprintf('Max element in first row = %f \n', max(A(1,:))));
disp(sprintf('Max element in first column = %f \n', max(A(:,1))));
if N > 20
    for tempCount = 1:20
        disp(sprintf('%e ', A(1, tempCount)));
    end
    disp(sprintf('\n'));
    for tempCount = 1:20
        disp(sprintf('%e ', A(tempCount, 1)));
    end
else
    disp(sprintf('%e ', A(1, tempCount)));
    disp(sprintf('%e ', A(tempCount, 1)));
end

%% Set up variables in Metropolis algorithm
movingParticleArr(totalSteps)=0;
pass = 0;
count_each_part(N)=0;
acceptCount = 0;

Nr_itr(totalSteps) = 0; % nr_in(totalSteps) = 0;
test_result(totalSteps) = 0; Bias_sign(totalSteps) = 0;
error_ratio(totalSteps) = 0;
reord = 0;
% match_global = zeros(size(A));
% Match_num = zeros(2, 1);
% Array_cut(1) = global_cutoff;
Array_dsum(1) = dsum;
% Match_Localnum = []; Match_num = [];
% Ssvl = [];

% kkk = 1;
%% start with same random seed
rng(0);
tic
for mcStep = 1:totalSteps
    % Generate a random move based upon normal distribution.
    [movingParticle, move] = random_move (N, sigma);
    randNum = rand;
    movingParticleArr(mcStep) = movingParticle;
    
    [u, ParticleConfig1] = move_particle (A, movingParticle, move, ParticleConfig, OrbitalConfig);
    
    fprintf(datafile, 'At step %d \n \tmoving part %d from %s to %s \n', mcStep,...
        movingParticle, pos_to_string(ParticleConfig(movingParticle,:)), pos_to_string(ParticleConfig1(movingParticle,:)));
    timeDense = tic;
    %--------------  Dense Method  --------------------
    detRatio_dense = dense_detRatio (A,u,movingParticle);
    timeDElapsed(mcStep) = toc(timeDense);
    %--------------  Iterative Method  -------------------
    %     [detRatio_iter, iter, err, flag, r_nrm] = fgmres_reorder (A, u, movingParticle, 1e-5);
    %     [detRatio_iter, iter, err, ~, ~, ~, ~, ~, in_iter] = fgmres_reorder_local(1e-3, A, u, movingParticle, 1e-4);
    timemcStep = tic;
    [detRatio_iter, mhat, iter, err] = iterative_detRatio(A, u, movingParticle, mcStep);
    detRatio = detRatio_dense;
    if (abs(err) > 6e-3)
        love = 1;
        Btemp = A;
    end
    if (detRatio*detRatio > randNum)
        fprintf(datafile, '-----accept   %f   %f \n\n', detRatio*detRatio, randNum);
        A(movingParticle,:) = A(movingParticle,:) + u';
        ParticleConfig = ParticleConfig1;
        acceptCount = acceptCount+1;
        globalCount = globalCount + 1;
        %% Using the variables for updating the preconditioner.
        array_m(:, end + 1) = mhat;
        array_u(:, end + 1) = u;
        array_detRatio(end + 1) = detRatio;
        %% reordering step
        [ A, dsum, global_cutoff, localcutoff, reord, Array_dsum, ParticleConfig] = combined_reorder(A, ...
            movingParticle, global_cutoff, ParticleConfig, reord, localcutoff, dsum, mcStep, Array_dsum);
    else
        fprintf(datafile, '-----reject    %f   %f \n\n', detRatio*detRatio, randNum);
    end
    timeElapsed = toc(timemcStep);
    timeIn(mcStep) = timeElapsed;
    
    Nr_itr(mcStep) = iter;
    %     nr_in(mcStep) = in_iter;
    error_ratio(mcStep) = err;
    %     f1 = detRatio_dense*detRatio_dense;
    %     f2 = detRatio_iter*detRatio_iter;
    f_sign = min(detRatio_dense*detRatio_dense, 1) - min(detRatio_iter*detRatio_iter, 1);
    %     f = abs(f_sign);
    test_result(mcStep) = abs(f_sign); Bias_sign(mcStep) = f_sign;
    if (abs(err) > 3e-3)
        fprintf(datafile1, 'At step %d, error is %8.4e and f is %8.4e  \n', mcStep, err, f_sign);
        if abs(f_sign) > 5e-3
            fprintf(datafile1, '\t f is larger than 5e-3!!!  \n');
        else
        end
    end
    if (abs(err) > 6e-3)
        temp1 = norm( Btemp - mat_L*mat_U, 'fro');
        temp2 = norm( eye(N) - mat_U\(mat_L\Btemp), 'fro');
        fprintf('At step %d, ||A-LU||_F is %8.4e and ||I-(LU)\A||_F is %8.4e  \n', mcStep, temp1, temp2);
        fprintf(datafile1, 'At step %d, ||A-LU||_F is %8.4e and ||I-(LU)\A||_F is %8.4e  \n', mcStep, temp1, temp2);
    else
    end
    fprintf(datafile, '\tdense: %f   iterative:%f \n', detRatio_dense, detRatio_iter);
    
    % Find if a pass is completed
    for count = 1:N
        count_each_part(count) = length(find(movingParticleArr == count));
    end
    % Check the observables every pass exactly
    if (rem(mcStep,N)== 0 && min(count_each_part) > pass)
        pass = pass + 1;
    end
    % Save configuration
    if (mod(mcStep, 5000) == 0)
        save(sprintf('%ddata', mcStep));
        fprintf(datafile1, sprintf('\n At step %d, accept move is %d, \n Accept ratio is %6.4f%%. \n', ...
            mcStep, acceptCount, 100*acceptCount/mcStep));
    end
end
timeAll = toc;
timeAllin = sum(timeIn);
timeDen = sum(timeDElapsed);
%%
fprintf(datafile1, '\n %d out of %d movements are accepted. Accept ratio is %6.4f%%  \n', ...
    acceptCount, mcStep, 100*acceptCount/mcStep);
fprintf(datafile1, '\n Finally, after %d steps, pass = %d. \n', totalSteps, pass);
fclose(datafile1);
fclose(datafile);
save('Finaldata');
