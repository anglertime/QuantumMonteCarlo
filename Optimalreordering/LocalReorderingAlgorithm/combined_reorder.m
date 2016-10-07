function [ A, globalCount, dsum, global_cutoff, localcutoff, reord, Array_cut, Array_dsum, ParticleConfig,...
    Match_num, Match_Localnum  ] = combined_reorder( A, movingParticle, global_cutoff, ParticleConfig, reord,...
    localcutoff, dsum, mcStep, globalCount, Match_Localnum, Array_cut, Array_dsum, Match_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global datafile1 reinitialize indicator

%% Initialization
N = size(A, 1);
if N < 500
    cycle = N;
else
    cycle = floor(0.6*N);
end
dsumcut = 1.2;
tempdsum = sum(diag(A));
accMovPdia = A(movingParticle, movingParticle);

%% Global reordering
if (tempdsum < dsum - dsumcut) || (mod(globalCount, cycle) == 0)
    if (tempdsum < dsum - dsumcut)
        fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 2. \n', mcStep));
    else
        fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 3. \n', mcStep));
    end
    globalCount = 1;    % Reorder globally
    indicator = ones(1, N);    % Initialize the indicator
    reinitialize = 1;
    reord = reord + 1;
%     Ssvl(reord) = min(svd(A));
    
    [ global_cutofftemp, ~] = max_min_diagonal( A, 2*N+2, 5e-3, 1e-3 );
    global_cutoff = min(global_cutoff, global_cutofftemp);
    [ A, ParticleConfig, ~, permu ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
    
    localcutoff = global_cutoff;
    Match_num(:, end + 1) = [mcStep; permu];
    dsum = sum(diag(A));  % Update dsum
    Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
    return;
else
end
%% Local reordering
if indicator(movingParticle) == 1
    if (accMovPdia < global_cutoff) % (tempdmin < min(0.8*cutoff, 8e-5))
        reinitialize = 1;   % Reorder locally
        reord = reord + 1;  
%         Ssvl(reord) = min(svd(A));
        
        [ A, ParticleConfig, localMatching, Templocalcutoff, localSize, Match_Localnum ] = LocalReorder( A, ParticleConfig, mcStep, Match_Localnum );
        dsum = sum(diag(A));   % Update dsum
        localcutoff = min(localcutoff, Templocalcutoff);
        Array_cut(reord + 1) = localcutoff; Array_dsum(reord + 1) = dsum;
        fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> global cut. Localsize %d \n', mcStep, movingParticle, localSize));
    else
    end
elseif (size(Match_Localnum, 1) > 1) && (accMovPdia < localcutoff) %&& (accMovPdia < 0.9*localcutoff)
    reinitialize = 1;   % Reorder locally
    reord = reord + 1;  
%     Ssvl(reord) = min(svd(A));
    
    [ A, ParticleConfig, localMatching, Templocalcutoff, localSize, Match_Localnum ] = LocalReorder( A, ParticleConfig, mcStep, Match_Localnum );
    dsum = sum(diag(A));   % Update dsum
    localcutoff = min(localcutoff, Templocalcutoff);
    Array_cut(reord + 1) = localcutoff; Array_dsum(reord + 1) = sum(diag(A));
    fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> local cut. Localsize %d \n', mcStep, movingParticle, localSize));
end

if (min(indicator) == -2)
    fprintf(datafile1, sprintf('\n At step %d, reordered the matrix globally because some particle has been permutated three times since one global reordering. \n', mcStep));
    globalCount = 1;    % Reorder globally
    indicator = ones(1, N);    % Initialize the indicator
    reinitialize = 1;
    reord = reord + 1;  
%     Ssvl(reord) = min(svd(A));
    
    [ global_cutofftemp, ~] = max_min_diagonal( A, 2*N+2, 5e-3, 1e-3 );
    global_cutoff = min(global_cutoff, global_cutofftemp);
    [ A, ParticleConfig, ~, permu ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
    
    localcutoff = global_cutoff;
    Match_num(:, end + 1) = [mcStep; permu];
    dsum = sum(diag(A));  % Update dsum
    Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
else
end
end