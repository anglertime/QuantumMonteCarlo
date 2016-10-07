function totalEnergy = Energy(coords, orbitals, M, gaussian_k, m, plank_const) 
% Find the total energy
laplaceTerm=0;
system_size = length(coords);
trans_inv_M = transpose(inv(M));
% For particles
for i=1:system_size
    % For orbitals
    for j=1:system_size 
        diff=coords(i,:)-orbitals(j,:);
        diff=put_in_box(diff); 
        dist=dot(diff,diff);
        laplaceTerm = laplaceTerm + (6*gaussian_k-4*(gaussian_k^2)*dist)*M(i,j)*trans_inv_M(i,j); 
    end
end
totalEnergy = plank_const^2 / (2*m*system_size) * laplaceTerm;