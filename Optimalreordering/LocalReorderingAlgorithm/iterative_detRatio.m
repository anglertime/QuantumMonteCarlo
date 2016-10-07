function  [detRatio, mhat, Iterations, err] = iterative_detRatio(A, u, movingParticle, mcStep, acceptCount)
% function  detRatio = iterative_detRatio(A, u, movingParticle)
% ILU-preconditioner and GMRES

global datafile detRatio_dense array_m array_u array_detRatio reinitialize mat_L mat_U
n = size(A, 2);
b = zeros(n, 1);
b(movingParticle) = 1;
tol = 5e-4;
maxit = n;

% Perform preconditioning every certain steps.
% if 1%( mcStep == 1 || rem(mcStep, 25) == 0 ) acceptCount
if ( mcStep == 1 || rem(acceptCount, 50) == 0 || reinitialize == 1)
    %         setup.type = 'nofill';
    %         [mat_L,mat_U] = ilu(sparse(A),setup);
    mat_L = []; mat_U = [];
    setup.type = 'ilutp';
%     setup.milu = 'row';
%     setup.type = 'crout';
    setup.droptol = 0.05;
    setup.thresh = 0.001;
    % Drop tolerance of the incomplete LU factorization. droptol is a non-negative scalar. The default value is 0, which produces the complete LU factorization.
    % The nonzero entries of U satisfy
    %   abs(U(i,j)) >= droptol*norm((A:,j)),
    % with the exception of the diagonal entries, which are retained regardless of satisfying the criterion. The entries of L are tested against the local drop tolerance before being scaled by the pivot, so for nonzeros in L
    % abs(L(i,j)) >= droptol*norm(A(:,j))/U(j,j).
    [mat_L,mat_U] = ilu(sparse(A),setup);
%     mat_L = P\mat_L;
    % Reinitialize these arrays
    array_m = [];
    array_u = [];
    array_detRatio = [];
    reinitialize = 0;
end

x0 = zeros(n, 1);
x0(movingParticle) = 1;

[mhat,flag,~,iter,resvec,~] = my_gmres(A,b,[],tol,maxit,mat_L,mat_U,x0);
% [mhat,flag,~,iter,resvec,~] = gmres(afun,b,[],tol,maxit);
% [mhat,flag,~,iter,resvec] = gmres(A,b,[],tol,maxit,@(x)mfun( x, A, P),[], x0);
Iterations = iter(1)*iter(2); % = length(resvec)-1.
% 
 fprintf(datafile, 'flag = %d, iter = %d, nnz(L) = %d, nnz(U) = %d, nnz(L+U) = %d, nnz(A) = %d \n',...
     flag, length(resvec)-1, nnz(mat_L), nnz(mat_U), nnz(mat_L+mat_U), nnz(A));
detRatio = 1+mhat'*u;
% err = abs(detRatio-detRatio_dense);
err = (detRatio-detRatio_dense);

function y = mfun( r, A, P)
%UNTITLED7 Summary of this function goes here
% two steps to express the prec. making use of matrix-vector product
% 1. v = P*r
% 2. v = U\(L\v)
y = P*r;
y = A*(U\(L\y));
end
end