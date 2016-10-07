function [B] = make_diagonal_dominant(B)

N = size(B,1);

for i=1:N
    
    % Get the index of max absolute value from the columns not yet touched.
    % That is from the later columns and not the before column.
    [maxValue, index] = max(abs(B(i,i:N)));
    index = index + i - 1;
    
    % If the index is the diagonal element, then check for the rows below
    % as we can move the particles too.
    if (index == i)
        [maxValue, index] = max(abs(B(i:N,i)));
        index = index + i - 1;
        
        if (index == i)
            % Do nothing.
        else
            % Swap the value in B matrix for those two rows.
            for innerCount = 1:N
                innerValue = B(index, innerCount);
                B(index, innerCount) = B(i, innerCount);
                B(i, innerCount) = innerValue;
            end
            
            % Swap the particle positions
%             tempArr = particleConfig1(i,:);
%             particleConfig1(i,:) = particleConfig1(index,:);
%             particleConfig1(index,:) = tempArr;
        end
    else
        % Swap the value in B matrix for those two columns.
        for innerCount = 1:N
            innerValue = B(innerCount, index);
            B(innerCount, index) = B(innerCount, i);
            B(innerCount, i) = innerValue;
        end

        % Swap the orbital position.
%         tempArr = orbitalConfig(i,:);
%         orbitalConfig(i,:) = orbitalConfig(index,:);
%         orbitalConfig(index,:) = tempArr;
    end
end

% Update u (assuming that at least one swap happened).
% for j=1:N
%      u(j, 1) = calculate_coef (particleConfig1(movingParticle,:), orbitalConfig(j,:)) - B(movingParticle, j);
% end
