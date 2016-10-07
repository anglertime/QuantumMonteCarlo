function [ A, dsum, global_cutoff, localcutoff, reord, Array_dsum, ParticleConfig] = combined_reorder( A, movingParticle, global_cutoff, ParticleConfig, reord,...
    localcutoff, dsum, mcStep, Array_dsum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global datafile1 reinitialize indicator globalCount

%% Initialization
N = size(A, 1);
% if N < 500
%     cycle = N;
% else
%     cycle = floor(0.6*N);
% end
if N < 2001
    dsumcut = 0.015*dsum;
else
    dsumcut = 0.03*dsum;
end
tempdsum = sum(diag(A));
accMovPdia = A(movingParticle, movingParticle);

%% Global reordering
if (tempdsum < dsum - dsumcut) %|| (mod(globalCount, cycle) == 0)
    globalCount = 1;    % Reorder globally
    indicator = ones(1, N);    % Initialize the indicator
    reinitialize = 1;
    reord = reord + 1;
%     Ssvl(reord) = min(svd(A));
    
    [ global_cutofftemp, ~] = max_min_diagonal( A, 2*N+2, 5e-3, 1e-3 );
    global_cutoff = min(global_cutoff, global_cutofftemp);
    [ A, ParticleConfig, ~, permu ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
    
    fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 2. Permutations %d. \n', mcStep, permu));
    localcutoff = global_cutoff;
    dsum = sum(diag(A));  % Update dsum
    Array_dsum(reord + 1) = dsum;
    return;
else
end
%% Local reordering
if indicator(movingParticle) == 1
    if (accMovPdia < 0.2*global_cutoff) % (tempdmin < min(0.8*cutoff, 8e-5))
        reord = reord + 1;  
        [ A, ParticleConfig, permu, Templocalcutoff, localSize] = LocalReorder( A, ParticleConfig);
        localcutoff = min(localcutoff, Templocalcutoff);
        fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> global cut. Localsize %d. Permutations %d. \n', mcStep, movingParticle, localSize, permu));
    else
    end
elseif (accMovPdia < 0.2*localcutoff) %&& (accMovPdia < 0.9*localcutoff)
    reord = reord + 1;  
    [ A, ParticleConfig, permu, Templocalcutoff, localSize] = LocalReorder( A, ParticleConfig);
    dsum = sum(diag(A));   % Update dsum
    localcutoff = min(localcutoff, Templocalcutoff);
    Array_dsum(reord + 1) = sum(diag(A));
    fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> local cut. Localsize %d. Permutations %d. \n', mcStep, movingParticle, localSize, permu));
end

if (min(indicator) == -3)
    fprintf(datafile1, sprintf('\n At step %d, reordered the matrix globally because some particle has been permutated four times since one global reordering. \n', mcStep));
    globalCount = 1;    % Reorder globally
    indicator = ones(1, N);    % Initialize the indicator
    reinitialize = 1;
    reord = reord + 1;  
    [ global_cutofftemp, ~] = max_min_diagonal( A, 2*N+2, 5e-3, 1e-3 );
    global_cutoff = min(global_cutoff, global_cutofftemp);
    [ A, ParticleConfig, ~, ~ ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
    
    localcutoff = global_cutoff;
    dsum = sum(diag(A));  % Update dsum
    Array_dsum(reord + 1) = dsum;
else
end
end