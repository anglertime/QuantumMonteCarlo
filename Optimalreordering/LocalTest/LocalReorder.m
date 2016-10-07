function [ B_re, newParticle, LocalMatching, cutoff, localSize, Match_LocalnumNew ] = LocalReorder( B, Particle, mcStep, Match_Localnum )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global max_Acut indicator
N = size(B, 1);
%% Pick out the local domain to reorder
% Make sure there is only one smallest positive element in the diagonal
% if length(find(diag(B) < 1e-5)) > 1
%     fprintf(datafile1, sprintf('\n At step %d, More than one element of the diagonal is smaller than 1e-5. \n', mcStep));
%     fprintf(datafile1, sprintf('\t Reorder the matrix globally. Flag = 5$$$$$$$$!!!!! \n'));
%     [ cutoff, ~] = max_min_diagonal( B, 2*N+2, 1e-3, 1e-3 );
%     [ B_re, newParticle, LocalMatching ] = max_weighted_match( B, cutoff, 3, Particle );
% %     Match_num(:, end + 1) = [mcStep; length(find(matching == 1:N))];
%     indicator = ones(1, N);    % Initialize the indicator
%     num = -1;
%     localSize = -1;
%     Match_LocalnumNew = Match_Localnum;
%     return;
% else
%     [~, CenterP] = min(diag(B));
% end
% Local submatrix has less choices so that the resulting mindiag
% should be smaller. However, we cant let it be too small!!!
% Orb = find(B(CenterP, :) >= B(CenterP, CenterP));
% dropTol = 1e-5;

%% Method 1.
[~, CenterP] = min(diag(B));
dropTol = 8*1/max_Acut;
Orb = find(B(CenterP, :) >= dropTol);
% Method 2
% Pick [CenterP - m, CenterP + m]
% width = 6;
% if (CenterP > width) && (CenterP < N - width)
%     Orb2 = (CenterP - width): (CenterP + width);
% elseif (CenterP <= width)
%     Orb2 = [1:CenterP, N - (width - CenterP) : N];
% elseif (CenterP >= N - width)
%     Orb2 = [1:N - CenterP, CenterP:N];
% end
% Orb = union(Orb, Orb2);

Par = Orb;
m1 = length(Par);
mapper = 1:N;
mapper(Par) = [];
mapper = [Par, mapper];
M = B(mapper, mapper);
LocalB = M(1:m1, 1:m1);

%%
[ cutoff, ~] = max_min_diagonal( LocalB, 2*m1+2, 1e-3, 1e-3 );
multiple = 3;% Scalar we use to replace elements less than the threshold
%Change the elements less than the threshold into larger ones so that we
% could discard them in the Hungarian algorithm
C = LocalB; replace = max(norm(LocalB, 1), norm(LocalB, inf));
C(C >= cutoff) = replace - C(C >= cutoff);
C(C < cutoff) = multiple*replace;

%% Use Hungarian Algorithm to find the maximum matching for B. We plug in C here
% because our Hungarian is to find the minimum weighted matching.
% [Matching, ~] = Hungarian(C);
% [part, orb] = find(Matching > 0); % Extract particle and orbital matching information
[LocalMatching, ~] = Munkres_Hungarian(C);

%% Reorder the matrix with obtained matching information
% permu = length(find( LocalMatching ~= 1:m1 ));
% % 2014-02-21 Ming. Create the mapping for input matrix
% mapping = 1:size(B,1);
% mapping(Matching(Matching ~= 1:size(B, 1))) = find( Matching ~= 1:size(B, 1) );
% if (nargin >= 4)
%     % Use Matching to remap the input matrix B
%     B_opt = B(mapping, :);
%     newParticle = particle(mapping, :);
% else
%     B_opt = B(mapping, :);
%     newParticle = [];
% end

%% Go through the matching
B_re = B;
permu = 0;
temp_permu = [];
for i = 1:m1
    if LocalMatching(i) == i
    else
        % For row Par(LocalMatching(i)), swap rows of row Par(LocalMatching(i)) and row Par(i).
        temp = B_re(Par(LocalMatching(i)), :);
        B_re(Par(LocalMatching(i)), :) = B_re(Par(i), :);
        B_re(Par(i), :) = temp;
        % Swap the particle positions
        tempArr = Particle(Par(LocalMatching(i)), :);
        Particle(Par(LocalMatching(i)), :) = Particle(Par(i), :);
        Particle(Par(i), :) = tempArr;
        % Swap particles orb(i) and part(i) in the stored matching
        Par(Par == Par(LocalMatching(i))) = Par(i);
        permu = permu + 1;
        temp_permu = union(temp_permu, [Par(LocalMatching(i)), Par(i)]);
    end
end
indicator(temp_permu) = indicator(temp_permu) - 1;    % Initialize the indicator for permutated rows
newParticle = Particle;
localSize = m1;
Match_LocalnumNew = Match_Localnum;
Match_LocalnumNew(1:3, end + 1) = [mcStep; permu; m1];
end

