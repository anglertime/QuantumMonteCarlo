function [ B_opt, cost, newParticle, Matching ] = max_weighted_matching( B, cutoff, particle )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
C = B;
C = (-min(200, size(B, 1))*(C < cutoff))+((C >= cutoff).*B);

[Matching, cost] = Hungarian(-C);

[part, orb] = find(Matching > 0);
B_opt = B;
for i = 1:size(B_opt, 1)
    % go through the matching
    % For orbital orb(i), swap part(i) into the diagonal. 
%     for innerCount = 1:n
%     end
    temp = B_opt(orb(i), :);
    B_opt(orb(i), :) = B_opt(part(i), :);
    B_opt(part(i), :) = temp;
     % Swap the particle positions
    tempArr = particle(orb(i),:);
    particle(orb(i),:) = particle(part(i),:);
    particle(part(i),:) = tempArr;
    
    % Swap particles orb(i) and part(i)
%     temp = part(part == orb(i));
    part(part == orb(i)) = part(i);
end
newParticle = particle;
end

