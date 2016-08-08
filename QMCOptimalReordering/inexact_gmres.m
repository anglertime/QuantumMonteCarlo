function [x,r_nrm,innertol,in_iter,iter,flag,in_all] = inexact_gmres(A11, A12, A21, A22, b,x0,tol,max_it,epsilon)

%  global A11 A12 A21 A22
 

   n = length(x0);             % initialization
   m = size(A22, 1);
   iter = 0;                                        
   flag = 0;
   r_nrm = zeros(max_it+1,1);

   x = x0;
   r0 = norm(A21*x0);   % This is when x0 = A11\e1.
   innertol(1) = epsilon/r0;
   iter = 1;
         %% Jacobian preconditioner
   dd = diag(A22);
   JcA22 = A22;
   for i = 1:m
       if dd(i) < 1e-2
           dd(i) = 1;
       else
           JcA22(i, :) = JcA22(i, :)/dd(i);
       end
   end
   %%
   Ax = mat_vec(x);
   r = b - Ax;
   r_nrm(1) = norm( r );

   if ( r_nrm(1) <= tol ) 
       in_all = sum(in_iter);
       flag = 1;
       return;
   end
   notconv = 1;

   %%
   i = 1;   % initialize outer steps
   V(1:n,1:i+1) = zeros(n,i+1);
   H(1:i+1,1:i) = zeros(i+1,i);
   cs(1:i) = zeros(i,1);
   sn(1:i) = zeros(i,1);
   e1    = zeros(i+1,1);
   e1(1) = 1.0;
   
   V(:,1) = r / r_nrm(iter);
   s = r_nrm(iter)*e1;       % initialize rhs Hessberg system
   
   iter = 2;

      while notconv && (iter - 1 <= max_it)
        w = mat_vec(V(:,i));
        for k = 1:i,
          H(k,i)= V(:,k)'*w;
          w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = norm( w );
        if H(i+1,i) > 0,
          V(:,i+1) = w / H(i+1,i);
          for kk = 1:i-1,                              % apply Givens rotation
            temp     =  cs(kk)*H(kk,i) + conj(sn(kk))*H(kk+1,i);
            H(kk+1,i) = -sn(kk)*H(kk,i) + conj(cs(kk))*H(kk+1,i);
            H(kk,i)   = temp;
          end
          [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
          temp   = cs(i)*s(i);                        % approximate residual norm
          s(i+1) = -sn(i)*s(i);
          s(i)   = temp;
          H(i,i) = cs(i)*H(i,i) + conj(sn(i))*H(i+1,i);
          H(i+1,i) = 0.0;
          r_nrm(iter)  = abs(s(i+1));

          if ( r_nrm(iter) <= tol ),
            notconv=0;
            flag = 0;
          else
            i = i+1;
          end % if ( r_nrm(iter) <= tol )
        else
           notconv=0;
        end % if H(i+1,i) > 0
        iter = iter+1;
      end % while notconv & (i<=m) & (iter <= max_it)

      
      if notconv == 1   %2013-04-10
          i = i - 1;
      end
      
      if iter == max_it
          flag = 2;
      end
      
      y = H(1:i,1:i) \ s(1:i);
      x = x + V(:,1:i)*y;
      
   iter = iter - 1;
   r_nrm = r_nrm(1:iter); % set to actual size
   in_all = sum(in_iter);
%%-------------------------------------------
    function y = mat_vec(x)
        %%     inner_tol = innertol(1);
        if iter == 1
            inner_tol = 1e-6;
            innertol(iter) = inner_tol;
        elseif (iter > 1) && (iter < 4)
            inner_tol = 1e-4/r_nrm(iter-1);
%             inner_tol = 1e-6/r_nrm(iter-1);
            %         inner_tol = 1e-4*(1/10)/r_nrm(iter-1);
            innertol(iter)= inner_tol;
        elseif (iter >= 4) %&& (iter < 6)
%             inner_tol = 1e-5/r_nrm(iter-1);
            inner_tol = 1e-4/r_nrm(iter-1);
            innertol(iter)= inner_tol;
%         else
%             inner_tol = 1e-4/r_nrm(iter-1);
%             innertol(iter)= inner_tol; 
        end
        z = A21*x;
        [sol,~,inner_iter,~] = mygmres(JcA22,z./dd,zeros(size(z)),inner_tol,size(A22, 1));
        inner_iter = inner_iter(1);
        in_iter(iter) = inner_iter;
        sol = A12*sol;
        y = A11*x - sol;
    end
end