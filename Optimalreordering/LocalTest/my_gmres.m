function [x,flag,relres,iter,resvec,reorder] = my_gmres(A,b,restart,tol,maxit,M1,M2,x0,varargin)
%GMRES   Generalized Minimum Residual Method.
%   X = GMRES(A,B) attempts to solve the system of linear equations A*X = B
%   for X.  The N-by-N coefficient matrix A must be square and the right
%   hand side column vector B must have length N. This uses the unrestarted
%   method with MIN(N,10) total iterations.
%
%   X = GMRES(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X) accepts a vector input X and returns the matrix-vector
%   product A*X. In all of the following syntaxes, you can replace A by
%   AFUN.
%
%   X = GMRES(A,B,RESTART) restarts the method every RESTART iterations.
%   If RESTART is N or [] then GMRES uses the unrestarted method as above.
%
%   X = GMRES(A,B,RESTART,TOL) specifies the tolerance of the method.  If
%   TOL is [] then GMRES uses the default, 1e-6.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT) specifies the maximum number of outer
%   iterations. Note: the total number of iterations is RESTART*MAXIT. If
%   MAXIT is [] then GMRES uses the default, MIN(N/RESTART,10). If RESTART
%   is N or [] then the total number of iterations is MAXIT.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M) and
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2) use preconditioner M or M=M1*M2
%   and effectively solve the system inv(M)*A*X = inv(M)*B for X. If M is
%   [] then a preconditioner is not applied.  M may be a function handle
%   returning M\X.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2,X0) specifies the first initial
%   guess. If X0 is [] then GMRES uses the default, an all zero vector.
%
%   [X,FLAG] = GMRES(A,B,...) also returns a convergence FLAG:
%    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
%    1 GMRES iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 GMRES stagnated (two consecutive iterates were the same).
%
%   [X,FLAG,RELRES] = GMRES(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL. Note with
%   preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   [X,FLAG,RELRES,ITER] = GMRES(A,B,...) also returns both the outer and
%   inner iteration numbers at which X was computed: 0 <= ITER(1) <= MAXIT
%   and 0 <= ITER(2) <= RESTART.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = GMRES(A,B,...) also returns a vector of
%   the residual norms at each inner iteration, including NORM(B-A*X0).
%   Note with preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = gmres(A,b,10,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to GMRES:
%      x1 = gmres(@(x)afun(x,n),b,10,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ, LUINC,
%   FUNCTION_HANDLE.

%   References
%   H.F. Walker, "Implementation of the GMRES Method Using Householder
%   Transformations", SIAM J. Sci. Comp. Vol 9. No 1. January 1988.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.21.4.9 $ $Date: 2007/02/02 23:21:42 $
reorder = 0;
if (nargin < 2)
   error('MATLAB:gmres:NumInputs','Not enough input arguments.');
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
   % Check matrix and right hand side vector inputs have appropriate sizes
   [m,n] = size(A);
   if (m ~= n)
      error('MATLAB:gmres:SquareMatrix','Matrix must be square.');
   end
   if ~isequal(size(b),[m,1])
      error('MATLAB:gmres:VectorSize', '%s %d %s', ...
         'Right hand side must be a column vector of length', m,...
         'to match the coefficient matrix.');
   end
else
   m = size(b,1);
   n = m;
   if ~isvector(b) || (size(b,2) ~= 1) % if ~isvector(b,'column')
      error('MATLAB:gmres:Vector','Right hand side must be a column vector.');
   end
end

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(restart) || (restart == n)
   restarted = false;
else
   restarted = true;
end
if (nargin < 4) || isempty(tol)
   tol = 1e-6;
end
if (nargin < 5) || isempty(maxit)
   if restarted
      maxit = min(ceil(n/restart),10);
   else
      maxit = min(n,10);
   end
end

if restarted
   outer = maxit;
   if restart > n
     warning('MATLAB:gmres:tooManyInnerIts', 'RESTART is %d %s\n%s %d.', ...
       restart, 'but it should be bounded by SIZE(A,1).', ...
       '         Setting RESTART to', n);
     restart = n;
   end
   inner = restart;
else
   outer = 1;
   if maxit > n
     warning('MATLAB:gmres:tooManyInnerIts', 'MAXIT is %d %s\n%s %d.', ...
       maxit, 'but it should be bounded by SIZE(A,1).', ...
       '         Setting MAXIT to', n);
     maxit = n;
   end
   inner = maxit;
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                      % Norm of rhs vector, b
if (n2b == 0)                       % if    rhs vector is all zeros
   x = zeros(n,1);                  % then  solution is all zeros
   flag = 0;                        % a valid solution has been obtained
   relres = 0;                      % the relative residual is actually 0/0
   iter = [0 0];                    % no iterations need be performed
   resvec = 0;                      % resvec(1) = norm(b-A*x) = norm(0)
   if (nargout < 2)
      itermsg('gmres',tol,maxit,0,flag,iter,NaN);
   end
   return
end

if ((nargin >= 6) && ~isempty(M1))
   existM1 = 1;
   [m1type,m1fun,m1fcnstr] = iterchk(M1);
   if strcmp(m1type,'matrix')
      if ~isequal(size(M1),[m,m])
         error('MATLAB:gmres:PreConditioner1Size', '%s %d %s',...
            'Preconditioner must be a square matrix of size', m, ...
            'to match the problem size.');
      end
   end
else
   existM1 = 0;
   m1type = 'matrix';
end

if ((nargin >= 7) && ~isempty(M2))
   existM2 = 1;
   [m2type,m2fun,m2fcnstr] = iterchk(M2);
   if strcmp(m2type,'matrix')
      if ~isequal(size(M2),[m,m])
         error('MATLAB:gmres:PreConditioner2Size', '%s %d %s', ...
            'Preconditioner must be a square matrix of size', m, ...
            'to match the problem size.');
      end
   end
else
   existM2 = 0;
   m2type = 'matrix';
end

if ((nargin >= 8) && ~isempty(x0))
   if ~isequal(size(x0),[n,1])
      error('MATLAB:gmres:XoSize', '%s %d %s', ...
         'Initial guess must be a column vector of length', n, ...
         'to match the problem size.');
   end
else
   x0 = zeros(n,1);
end
x = x0;

if ((nargin > 8) && strcmp(atype,'matrix') && ...
      strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
   error('MATLAB:gmres:TooManyInputs','Too many input arguments.');
end

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                 % Norm of residual

if (normr <= tolb)               % Initial guess is a good enough solution
   flag = 0;
   relres = normr / n2b;
   iter = [0 0];
   resvec = normr;
   if (nargout < 2)
      itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
   end
   return
end

if existM1
   r1 = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
   if any(~isfinite(r1))
      flag = 2;
      x = xmin;
      relres = normr / n2b;
      iter = [0 0];
      resvec = normr;
      return
   end
else
   r1 = r;
end

if existM2
   r = iterapp('mldivide',m2fun,m2type,m2fcnstr,r1,varargin{:});
   if any(~isfinite(r))
      flag = 2;
      x = xmin;
      relres = normr / n2b;
      iter = [0 0];
      resvec = normr;
      return
   end
else
   r = r1;
end
% Modified GMRES with rank-one update - Kapil, July 13th 2008.
r = update_precond(r);

normr = norm(r);
resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = zeros(2,inner);

% For reorder: Kapil 9th Jan 2010
% reorder = 0;

for outiter = 1 : outer
   %  Construct u for Householder reflector.
   %  u = r + sign(r(1))*||r||*e1
   u = r; beta = scalarsign(r(1))*normr;
   u(1) = u(1) + beta;
   u = u / norm(u);

   U = zeros(n,inner);
   Psol = zeros(n,inner);
   R = zeros(n,inner);
   U(:,1) = u;

   %  Apply Householder projection to r.
   %  w = r - 2*u*u'*r;
   w = zeros(n,1); w(1) = -beta;

   for initer = 1 : inner
      %  Form P1*P2*P3...Pj*ej.
      %  v = Pj*ej = ej - 2*u*u'*ej
      v = -2*(u(initer)')*u;
      v(initer) = v(initer) + 1;
      for k = (initer-1):-1:1
         v = v - 2*U(:,k)*(U(:,k)'*v);
      end

      Psol(:,initer) = v;

      % For reorder: Kapil 9th Jan 2010
      % Hopefully Matlab creates copy of variable.
%       old_v = v;   

      %  Apply A to v.
      v = iterapp('mtimes',afun,atype,afcnstr,v,varargin{:});
      %  Apply Preconditioner.
      if existM1
         v1 = iterapp('mldivide',m1fun,m1type,m1fcnstr,v,varargin{:});
         if any(~isfinite(v1))
            flag = 2;
            break
         end
      else
         v1 = v;
      end

      if existM2
         v = iterapp('mldivide',m2fun,m2type,m2fcnstr,v1,varargin{:});
         if any(~isfinite(v))
            flag = 2;
            break
         end
      else
         v = v1;
      end
      % Modified GMRES with rank-one update - Kapil, July 13th 2008.
      v = update_precond(v);

      % For reorder: Kapil 9th Jan 2010
%       v_norm = norm (old_v - v) / norm (old_v);
%       if (v_norm >= reorder)
%         reorder = v_norm;
%       end
      
      %  Form Pj*Pj-1*...P1*Av.
      for k = 1:initer
         v = v - 2*U(:,k)*(U(:,k)'*v);
      end

      %  Determine Pj+1.
      if (initer ~= length(v))
         %  Construct u for Householder reflector Pj+1.
         u = zeros(n,1);
         vhat = v(initer+1:end);
         alpha = norm(vhat);
         if (alpha ~= 0)
             alpha = scalarsign(vhat(1))*alpha;
             %  u = v(initer+1:end) +
             %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
             u(initer+1:end) = vhat;
             u(initer+1) = u(initer+1) + alpha;
             u = u / norm(u);
             U(:,initer+1) = u;

             %  Apply Pj+1 to v.
             %  v = v - 2*u*(u'*v);
             v(initer+2:end) = 0;
             v(initer+1) = -alpha;
         end
      end

      %  Apply Given's rotations to the newly formed v.
      for colJ = 1:initer-1
         tmpv = v(colJ);
         v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
         v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
      end

      %  Compute Given's rotation Jm.
      if ~(initer==length(v))
         rho = norm(v(initer:initer+1));
         J(:,initer) = v(initer:initer+1)./rho;
         w(initer+1) = -J(2,initer).*w(initer);
         w(initer) = conj(J(1,initer)).*w(initer);
         v(initer) = rho;
         v(initer+1) = 0;
      end

      R(:,initer) = v;

      if initer < inner
         normr = abs(w(initer+1));
         resvec( (outiter-1)*inner+initer+1 ) = normr;
      end

      if normr <= normrmin
         normrmin = normr;
         imin = outiter;
         jmin = initer;
      end

      if normr < tolb
         flag = 0;
         iter = [outiter, initer];
         break
      end
   end         % ends innner loop

   y = R(1:jmin,:) \ w(1:jmin);
   additive = Psol*y;
   x = x + additive;
   xmin = x;

   r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});

   if existM1
      r1 = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
      if any(~isfinite(r1))
         flag = 2;
         break
      end
   else
      r1 = r;
   end

   if existM2
      r = iterapp('mldivide',m2fun,m2type,m2fcnstr,r1,varargin{:});
      if any(~isfinite(r))
         flag = 2;
         break
      end
   else
      r = r1;
   end
   % Modified GMRES with rank-one update - Kapil, July 13th 2008.
   r = update_precond(r);

   normr = norm(r);

   resvec((outiter-1)*inner+initer+1) = normr;

   if normr <= normrmin
      xmin = x;
      normrmin = normr;
   end

   %  Test for stagnation on outer iterate.
   if flag~=2
      stagtest = zeros(n,1);
      ind = (x ~=0 );
      stagtest(ind) = additive(ind) ./ x(ind);
      if ( norm(stagtest,inf) < eps )
          flag = 3;
          break;
         %  No change in outer iterate.
      end
   end

   if normr < tolb
      flag = 0;
      iter = [outiter, initer];
      break;
   end
end         % ends outer loop

% returned solution is that with minimum residual
if flag == 0
   relres = normr / n2b;
else
   x = xmin;
   iter = [imin jmin];
   relres = normrmin / n2b;
end

% truncate the zeros from resvec
if flag <= 1 || flag == 3
   resvec = resvec(1:(outiter-1)*inner+initer+1);
   indices = resvec==0;
   resvec = resvec(~indices);
else
   if initer == 0
      resvec = resvec(1:(outiter-1)*inner+1);
   else
      resvec = resvec(1:(outiter-1)*inner+initer);
   end
end

% only display a message if the output flag is not used
if nargout < 2
   if restarted
      itermsg(sprintf('gmres(%d)',restart),tol,maxit,[i j],flag,iter,relres);
   else
      itermsg(sprintf('gmres'),tol,maxit,j,flag,iter(2),relres);
   end
end


function sgn = scalarsign(d)
sgn = sign(d);
if (sgn == 0)
    sgn = 1;
end


% Modified GMRES with rank-one update - Kapil, July 13th 2008.
% Note: MATLAB passes data only by value, also array's are also 'Variables'
% in definition of Matlab.
function vec = update_precond(vec)
global array_m array_u array_detRatio

for count = 1:length(array_detRatio)
    vec = vec - ((array_u(:,count)'*vec)/array_detRatio(count)) * array_m(:,count);
end