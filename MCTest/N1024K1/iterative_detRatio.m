function  [detRatio, mhat, Iterations] = iterative_detRatio(A, u, movingParticle, mcStep)
% function  detRatio = iterative_detRatio(A, u, movingParticle)
% ILU-preconditioner and GMRES

global datafile array_m array_u array_detRatio 
global globalCount reinitialize mat_L mat_U
n = size(A, 2);
b = zeros(n, 1);
b(movingParticle) = 1;
tol = 6e-4;
maxit = 200;

% Perform preconditioning every certain steps.
% if 1%( mcStep == 1 || rem(mcStep, 25) == 0 ) acceptCount
if ( mcStep == 1 || rem(globalCount, 62) == 0 || reinitialize == 1)
    %         setup.type = 'nofill';
    %         [mat_L,mat_U] = ilu(sparse(A),setup);
    %     mat_L = []; mat_U = [];
    setup.type = 'ilutp';
    %     setup.milu = 'row';
    %     setup.type = 'crout';
    % 0.01 for smaller system
    setup.droptol = 0.005;
    setup.thresh = 0.005;
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
    globalCount = 1;
end

% x0 = zeros(n, 1);
% x0(movingParticle) = 1;
[mhat,flag,~,iter,~,~] = my_gmres(A,b,[],tol,maxit,mat_L,mat_U,b);
% [mhat,flag,~,iter,resvec,~] = gmres(afun,b,[],tol,maxit);
% [mhat,flag,~,iter,resvec] = gmres(A,b,[],tol,maxit,@(x)mfun( x, A, P),[], x0);
Iterations = iter(1)*iter(2); % = length(resvec)-1.
detRatio = 1+mhat'*u;
% err = abs(detRatio-detRatio_dense);
% err = (detRatio-detRatio_dense);
% if (abs(err) > 3e-3)
%     temp1 = norm( A - mat_L*mat_U, 'fro');
%     temp2 = norm( eye(1024) - mat_U\(mat_L\A), 'fro');
%     fprintf('At step %d, ||A-LU||_F is %8.4e and ||I-(LU)\A||_F is %8.4e  \n', mcStep, temp1, temp2);
%     fprintf(datafile1, 'At step %d, ||A-LU||_F is %8.4e and ||I-(LU)\A||_F is %8.4e  \n', mcStep, temp1, temp2);
% else
% end
% fprintf(datafile, 'flag = %d, iter = %d, nnz(L) = %d, nnz(U) = %d, nnz(L+U) = %d, nnz(A) = %d, err = %8.4e \n',...
%     flag, length(resvec)-1, nnz(mat_L), nnz(mat_U), nnz(mat_L+mat_U), nnz(A), err);
fprintf(datafile, 'flag = %d, iter = %d, nnz(L+U) = %d, nnz(A) = %d \n',...
    flag, Iterations, nnz(mat_L+mat_U), nnz(A));
end