function [ B_opt, newParticle, Matching, num ] = max_weighted_match_test( B, cutoff, multiple, particle, mcStep)
% Maximum weighted matching with threshold of diagonal elements
%
%   Input
%     B          :      input matrix. Notice that all elements should be positive.
%     cutoff     :      threshold of diagonal elments
%     particle   :      particle configuration
%
%   Output
%     B_opt      :      new matrix after we reorder the input matrix with maximum weighted matching
%     cost       :      cost of the Hungarian algorithm to find the matching
%     newParticle:      new particle configuration after reordering
%     Matching   :      the desired matching
%%
global datafile1
if (nargin < 2)
    error('Not enough inputs!');
elseif (nargin < 3)
    multiple = 3;% Scalar we use to replace elements less than the threshold
end
%Change the elements less than the threshold into larger ones so that we
% could discard them in the Hungarian algorithm
C = B; replace = max(norm(B, 1), norm(B, inf));
C(C >= cutoff) = replace - C(C >= cutoff);
C(C < cutoff) = multiple*replace;

% Use Hungarian Algorithm to find the maximum matching for B. We plug in C here
% because our Hungarian is to find the minimum weighted matching.
% [Matching, ~] = Hungarian(C);
% [part, orb] = find(Matching > 0); % Extract particle and orbital matching information
[Matching, ~] = Munkres_Hungarian(C);

%% Reorder the matrix with obtained matching information
% num = length(find( Matching ~= 1:size(B, 1) ));
% if (nargin == 4)
%     % Go through the matching
%     B_opt = B(Matching, :);
%     newParticle = particle(Matching, :);
% else
%     B_opt = B(Matching, :);
%     newParticle = [];
% end
%%
part = 1:size(B, 1); orb = Matching; num = 0;
B_opt = B;
if (nargin >= 4)
    % Go through the matching
    for i = 1:size(B_opt, 1)
        % For orbital orb(i), swap rows of part(i) and orb(i).
        if orb(i) == part(i) % part(i) is actually equal to i initially
        else
            fprintf(datafile1, sprintf('\n At step %d, swap row %d and row %d. \n', mcStep, orb(i), part(i)));
            temp = B_opt(orb(i), :);
            B_opt(orb(i), :) = B_opt(part(i), :);
            B_opt(part(i), :) = temp;
            % Swap the particle positions
            tempArr = particle(orb(i),:);
            particle(orb(i),:) = particle(part(i),:);
            particle(part(i),:) = tempArr;
            % Swap particles orb(i) and part(i) in the stored matching
            part(part == orb(i)) = part(i);
            num = num + 1;
        end
    end
    newParticle = particle;
else
    % Go through the matching
    for i = 1:size(B_opt, 1)
        if orb(i) == part(i)
        else
            % For orbital orb(i), swap rows of part(i) and orb(i).
            temp = B_opt(orb(i), :);
            B_opt(orb(i), :) = B_opt(part(i), :);
            B_opt(part(i), :) = temp;
            % Swap particles orb(i) and part(i) in the stored matching
            part(part == orb(i)) = part(i);
            num = num + 1;
        end
    end
    newParticle = [];
end
end

