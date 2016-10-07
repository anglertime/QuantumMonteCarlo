% clear global;
% clear;
% clear error_ratio
% clear Nr_itr
% load('Equi3DN1024K1MC15K.mat')
global constK Box datafile dim detRatio_dense sideLength max_dist_sqr max_Acut datafile1
global indicator reinitialize array_m  array_u  array_detRatio

%Initialization
dim = 3; N = 1024; indicator = ones(1, N); totalSteps = 3000;

% sigma = Box(1)*0.017; max_dist_sqr = max_dist;
[ParticleConfig, OrbitalConfig, Box, sigma, constK, sideLength, max_dist_sqr, max_Acut] = read_configuration (N);
% A = compute_matrix (ParticleConfig, OrbitalConfig);
load('30000A', 'A');
%% Initialization
[ global_cutoff, ~ ] = max_min_diagonal ( A, 2*N+2, 1e-3, 1e-3);
localcutoff = global_cutoff;

reinitialize = 0;
dsum = sum(diag(A)); dsumcut = 1; % 1.2 > 1.
if N < 500
    cycle = N;
else
    cycle = floor(0.8*N);
end

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

Nr_itr(totalSteps) = 0; nr_in(totalSteps) = 0;
test_result(totalSteps) = 0; Bias_sign(totalSteps) = 0;
error_ratio(totalSteps) = 0;

reord = 0; match_global = zeros(size(A));
% Match_num = zeros(2, 1);
Array_cutG(1) = global_cutoff; 
Array_cutL(1) = global_cutoff; Array_dsum(1) = dsum;
globalCount = 0; Match_Localnum = []; Match_num = [];
Ssvl = [];

kkk = 1;
%% start with same random seed
rng(0);
for mcStep = 1:totalSteps
    % Generate a random move based upon normal distribution.
    [movingParticle, move] = random_move (N, sigma);
    randNum = rand;
    movingParticleArr(mcStep) = movingParticle;
    
    [u, ParticleConfig1] = move_particle (A, movingParticle, move, ParticleConfig, OrbitalConfig);
    
    fprintf(datafile, 'At step %d \n \tmoving part %d from %s to %s \n', mcStep,...
        movingParticle, pos_to_string(ParticleConfig(movingParticle,:)), pos_to_string(ParticleConfig1(movingParticle,:)));
    
    %--------------  Dense Method  --------------------
    detRatio_dense = dense_detRatio (A,u,movingParticle);
    
    %--------------  Iterative Method  -------------------
    %     [detRatio_iter, iter, err, flag, r_nrm] = fgmres_reorder (A, u, movingParticle, 1e-5);
    %     [detRatio_iter, iter, err, ~, ~, ~, ~, ~, in_iter] = fgmres_reorder_local(1e-3, A, u, movingParticle, 1e-4);
    [detRatio_iter, mhat, iter, err] = iterative_detRatio(A, u, movingParticle, mcStep, acceptCount);
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
    fprintf(datafile, '\tdense: %f   iterative:%f \n', detRatio_dense, detRatio_iter);
    
    detRatio = detRatio_dense;
    if (detRatio*detRatio > randNum)
        fprintf(datafile, '-----accept   %f   %f \n\n', detRatio*detRatio, randNum);
        A(movingParticle,:) = A(movingParticle,:) + u';
        ParticleConfig = ParticleConfig1;
        acceptCount = acceptCount+1;
        globalCount = globalCount + 1;
        
        %% reordering step
        tempdmin = min(diag(A));
        accMovPdia = A(movingParticle, movingParticle);
        tempdsum = sum(diag(A));
        %         [ A, globalCount, dsum, global_cutoff, localcutoff, reord, Ssvl, Array_cut, Array_dsum, ParticleConfig,...
        %         Match_num, Match_Localnum ] = combined_reorder(A, movingParticle, global_cutoff, ParticleConfig, reord,...
        %         localcutoff, accMovPdia, tempdsum, dsum, dsumcut, cycle, mcStep, globalCount, Match_Localnum, Ssvl, Array_cut, Array_dsum, Match_num );
        
        %         if (tempdsum < dsum - dsumcut) || (mod(globalCount, cycle) == 0) || (tempdmin < 0.6*global_cutoff)
        %             [ global_cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
        %             [ A, ParticleConfig, matching ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
        %             dsum = sum(diag(A));  % Update dsum
        %             globalCount = 1;    % Reorder globally
        %             reinitialize = 1;
        %             reord = reord + 1;  Ssvl(reord) = min(svd(A));
        %             Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
        %             fprintf(datafile1, sprintf('\n At step %d, reordered the matrix globally. \n', mcStep));
        %         else
        %         end
        if (tempdmin < global_cutoff)
            [ global_cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
            [ B, ParticleConfig2, matching, swapNum ] = max_weighted_match_test( A, global_cutoff, 3, ParticleConfig, mcStep );
            reord = reord + 1;  % Reorder globally
            if reord <= N
                match_global(reord, :) = matching;
            end
            difnu = find(matching~=1:N);
            Match_num(reord) = swapNum;
            dsum = sum(diag(B));  % Update dsum
            globalCount = 1;    
            reinitialize = 1;
            Ssvl(reord) = min(svd(B));
            Array_cutG(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
            fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix globally. \n\n\n', mcStep, movingParticle));
            
            % Reorder locally
            [ C, ParticleConfig3, match_local, localMatching_num, Templocalcutoff, localSize, Match_Localnum, rowPerm ] = LocalReorder_chk( A, ParticleConfig, mcStep, Match_Localnum );
            dsum = sum(diag(C));   % Update dsum
            SamePermRownum(reord) = length(intersect(difnu, rowPerm));
            localcutoff = min(localcutoff, Templocalcutoff);
            Array_cutL(reord + 1) = localcutoff; Array_dsum(reord + 1) = dsum;
            fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> gglobal cut. Localsize %d \n\n\n', mcStep, movingParticle, localSize));
            if length(swapNum) <= 10
                if length(swapNum) >= 3
                break;
                end
            end
            A = B;
        else
        end
        
        if Array_cutL(end) <= 1/max_Acut
            fprintf(datafile1, sprintf('\n The cutoff is ZERO. \n'));
        end
        %% Using the variables for updating the preconditioner.
        array_m(:, end + 1) = mhat;
        array_u(:, end + 1) = u;
        array_detRatio(end + 1) = detRatio;
        %         if (mod(acceptCount, 100) == 0) && (acceptCount > 1)
        % %             cond1(kkk) = cond(A);
        % %             Ssvl(kkk) = svds(A, 1, 0);
        %             Ssvl(kkk) = min(svd(A));
        %             kkk = kkk + 1;
        %         end
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
    end
    % Save configuration
    if (mod(mcStep, 5000) == 0)
        %         Steps = fopen(sprintf('%dparticle.txt', mcStep), 'wt+');
        %         save_config(Steps, N, ParticleConfig);
        %         fclose(Steps);
        save(sprintf('%ddata', mcStep));
        fprintf(datafile1, sprintf('\n At step %d, accept move is %d, \n Accept ratio is %6.4f%%. \n', ...
            mcStep, acceptCount, 100*acceptCount/mcStep));
    end
end
%%
fprintf(datafile1, '\n %d out of %d movements are accepted. Accept ratio is %6.4f%%  \n', ...
    acceptCount, mcStep, 100*acceptCount/mcStep);
fprintf(datafile1, '\n Finally, after %d steps, pass = %d. \n', totalSteps, pass);
fclose(datafile1);
fclose(datafile);

save('Finaldata');
