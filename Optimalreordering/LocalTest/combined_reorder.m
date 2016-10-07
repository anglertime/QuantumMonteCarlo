function [ A, globalCount, dsum, global_cutoff, localcutoff, reord, Ssvl, Array_cut, Array_dsum, ParticleConfig,...
    Match_num, Match_Localnum  ] = combined_reorder( A, movingParticle, global_cutoff, ParticleConfig, reord,...
    localcutoff, accMovPdia, tempdsum, dsum, dsumcut, cycle, mcStep, globalCount, Match_Localnum, Ssvl, Array_cut, Array_dsum, Match_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global datafile1 reinitialize indicator
N = size(A, 1); cycle = 0.5*cycle;
if indicator(movingParticle) == 1
    if (accMovPdia < 0.9*global_cutoff) % (tempdmin < min(0.8*cutoff, 8e-5))
        if (tempdsum < dsum - dsumcut) || (mod(globalCount, cycle) == 0)
            if (tempdsum < dsum - dsumcut)
                fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 1a. \n', mcStep));
            else
                fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 1b. \n', mcStep));
            end
            globalCount = 1;    % Reorder globally
            indicator = ones(1, N);    % Initialize the indicator
            reinitialize = 1;
            reord = reord + 1;  Ssvl(reord) = min(svd(A));
            [ global_cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
            localcutoff = global_cutoff;
            [ A, ParticleConfig, matching ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
            Match_num(:, end + 1) = [mcStep; length(find(matching ~= 1:N))];
            dsum = sum(diag(A));  % Update dsum
            Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
        else
            reinitialize = 1;   % Reorder locally
            reord = reord + 1;  Ssvl(reord) = min(svd(A));
            [ A, ParticleConfig, localMatching, localMatching_num, Templocalcutoff, localSize, Match_Localnum ] = LocalReorder( A, ParticleConfig, mcStep, Match_Localnum );
            dsum = sum(diag(A));   % Update dsum
            localcutoff = min(localcutoff, Templocalcutoff);
            Array_cut(reord + 1) = localcutoff; Array_dsum(reord + 1) = dsum;
            fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> gglobal cut. Localsize %d \n', mcStep, movingParticle, localSize));
        end
    elseif (tempdsum < dsum - dsumcut) || (mod(globalCount, cycle) == 0)
        if (tempdsum < dsum - dsumcut)
            fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 2. \n', mcStep));
        else
            fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 3. \n', mcStep));
        end
        globalCount = 1;
        reinitialize = 1;   % Reorder globally
        indicator = ones(1, N);    % Initialize the indicator
        reord = reord + 1;  Ssvl(reord) = min(svd(A));
        [ global_cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
        localcutoff = global_cutoff;
        [ A, ParticleConfig, matching ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
        Match_num(:, end + 1) = [mcStep; length(find(matching ~= 1:N))];
        dsum = sum(diag(A));   % Update dsum
        Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
    end
elseif (size(Match_Localnum, 1) > 1) && (accMovPdia < 0.9*localcutoff) %&& (accMovPdia < 0.9*global_cutoff)
    if (tempdsum < dsum - dsumcut) || (mod(globalCount, cycle) == 0)
        if (tempdsum < dsum - dsumcut)
            fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 2loc. \n', mcStep));
        else
            fprintf(datafile1, sprintf('\n At step %d, reordered the matrix because flag = 3loc. \n', mcStep));
        end
        globalCount = 1;    % Reorder globally
        indicator = ones(1, N);    % Initialize the indicator
        reinitialize = 1;
        reord = reord + 1;  Ssvl(reord) = min(svd(A));
        [ global_cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
        localcutoff = global_cutoff;
        [ A, ParticleConfig, matching ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
        Match_num(:, end + 1) = [mcStep; length(find(matching ~= 1:N))];
        dsum = sum(diag(A));  % Update dsum
        Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
    else
        reinitialize = 1;   % Reorder locally
        reord = reord + 1;  Ssvl(reord) = min(svd(A));
        [ A, ParticleConfig, localMatching, localMatching_num, Templocalcutoff, localSize, Match_Localnum ] = LocalReorder( A, ParticleConfig, mcStep, Match_Localnum );
        dsum = sum(diag(A));   % Update dsum
        localcutoff = min(localcutoff, Templocalcutoff);
        Array_cut(reord + 1) = localcutoff; Array_dsum(reord + 1) = sum(diag(A));
        fprintf(datafile1, sprintf('\n At step %d, moveP %d, reordered the matrix LOCALLY -> local cut. Localsize %d \n', mcStep, movingParticle, localSize));
    end
end

if (min(indicator) == -2)
    fprintf(datafile1, sprintf('\n At step %d, reordered the matrix globally because some particle has been permutated three times since one global reordering. \n', mcStep));
    globalCount = 1;    % Reorder globally
    indicator = ones(1, N);    % Initialize the indicator
    reinitialize = 1;
    reord = reord + 1;  Ssvl(reord) = min(svd(A));
    [ global_cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
    localcutoff = global_cutoff;
    [ A, ParticleConfig, matching ] = max_weighted_match( A, global_cutoff, 3, ParticleConfig );
    Match_num(:, end + 1) = [mcStep; length(find(matching ~= 1:N))];
    dsum = sum(diag(A));  % Update dsum
    Array_cut(reord + 1) = global_cutoff; Array_dsum(reord + 1) = dsum;
else
end
end