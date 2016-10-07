t = 77; M = zeros(size(A)); for i = 1:N; M(i, Match(t, i)) = 1; end; spy(M, 'r*');
figure(2); spy(M - diag(diag(M)), 'b*');
min(diag(M))