function  [detRatio, iter, error, err1, local_flag, r_nrm, innertol, inner_nr, in_iter_all] = fgmres_reorder_local(drop_co, A, u, movingParticle, epsln)
% two level domain decomposition-preconditioner and GMRES

global datafile detRatio_dense

% Choose orbitals according to the realationship between the corresponding
% elements in matrix A. Local orbitals are the union of Ak and Ak+1.
B = A;
B(movingParticle,:) = B(movingParticle,:) + u';  % Ak+1 = B

[max_coe, max_index] = max(abs(A(movingParticle, :)));
[max_coeb, max_indexb] = max(abs(B(movingParticle, :)));
% Define dropping off coeficient
% drop_co = 1e-3;
V1a = find(abs(A(movingParticle, :) > max_coe*drop_co));
V1b = find(abs(B(movingParticle, :) > max_coeb*drop_co));
V1 = union(V1a, V1b);
% V1 = sort(V1);     % default is ascend
V1p = V1;
% This guarentee that the moving particle row is in the domain 1
if isempty(find(V1 == movingParticle, 1))
    % If the moving particle is not in domain 1, select the row/particle
    V1 = union(V1, movingParticle);
    V1p = V1;
    local_flag = 1;
    fprintf(datafile, '\n\n\t The moving index orbital is not include in the cutoff radius of moving particle\n');
else
    local_flag = 0;
end

m = length(V1p);
n = size(A, 2);

% Initial the matrix mapper
mapper_p = zeros(1,n); 
mapper_o = zeros(1,n); 
for i = 1:n
    mapper_p(i) = i;
    mapper_o(i) = i;
end
mapper_p(V1p) = [];
mapper_p = [V1p, mapper_p];
mapper_o(V1) = [];
mapper_o = [V1, mapper_o];

% Reordering
M = A(mapper_p, mapper_o);
A11 = M(1:m,1:m);
A12 = M(1:m,m+1:n);
A21 = M(m+1:n,1:m);
A22 = M(m+1:n,m+1:n);
% cond_A11 = cond(A11); cond_A22 = cond(A22);
for i = 1:m
    if V1p(i) == movingParticle
        pointer = i;
    end
end

% Get the fixed point solution
ek = zeros(m,1); ek(pointer) = 1;
u = u(mapper_o);
u1 = u(1:m, 1);
u2 = u(m+1:n, 1);
% x0 = zeros(size(ek));
x0 = A11\ek;
fprintf(datafile, '\tLocal domain is %d x %d  \n', size(A11, 1), size(A11, 2));

[z1,r_nrm,innertol,inner_nr,iter,flag,in_iter_all] = inexact_gmres(A11, A12, A21, A22, ek,x0,1e-6,length(x0),epsln);
% [z1,flag,r_nrm,iter,~] = gmres(@(x)mat_vec(x),ek,[],1e-6);
% iter = iter(2);
fprintf(datafile, 'flag = %d ; iter = %d ; inexact inner GMRES: %d, ||z1-x0||=%8.4e, ||u2|| is %8.4e. \n',...
    flag, iter, in_iter_all, norm(z1-x0), norm(u2));

detRatio = 1+u1'*z1;
detRatio1 = 1+u1'*x0;
error = detRatio - detRatio_dense;
% error = abs(error);
err1 = abs(detRatio1-detRatio_dense);
fprintf(datafile, 'Using x0=A11^-1ek: detRatio1=%8.4e   detRatio_err=%8.4e \n', detRatio1, err1);
fprintf(datafile, 'Inexact GMRES u1^T*z1:  detRatio=%8.4e   detRatio_err=%8.4e \n', detRatio, error);