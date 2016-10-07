Orb = [];
Orb2 = [];
temp1 = find(A(:, 271) >= 1e-5);
temp2 = find(A(271, :) >= 1e-5);
Orb = union(temp1, temp2);
for i = 1:length(Orb)
    temp = find(A(:, Orb(i)) > 1e-5);
    Orb2 = union(Orb2, temp);
end
length(intersect(Orb2, Orb))
length(intersect(Orb2, permu))
length(intersect(Orb, permu))
