function  [detRatio, iter, error, flag, r_nrm] = fgmres_reorder(A, u, movingParticle, tol)
% two level domain decomposition-preconditioner and GMRES

global datafile detRatio_dense

n = size(A, 1);
x0 = zeros(n, 1);
b= x0; b(movingParticle) = 1;

dd = diag(A);
JcA = A;
for i = 1:n
    if dd(i) < 1e-2
        dd(i) = 1;
    else
        JcA(i, :) = JcA(i, :)/dd(i);
    end
end

[z, r_nrm, iter, flag] = mygmres(JcA, b./dd, x0, tol, n); 

fprintf(datafile, 'flag = %d ; iter = %d \n', flag, iter);

detRatio = 1+u'*z;
error = detRatio - detRatio_dense;
error = abs(error);
fprintf(datafile, 'Inexact GMRES u1^T*z1:  detRatio=%8.4e   detRatio_err=%8.4e \n', detRatio, error);