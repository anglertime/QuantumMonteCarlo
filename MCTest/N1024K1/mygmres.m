function [x,r_nrm,iter,flag] = mygmres(mats,b,x0,tol,max_it)
   
   n = length(x0);             % initialization
   iter = 0;                                        
   flag = 0;
   r_nrm = zeros(max_it+1,1);

   x = x0;
   Ax = mats*x;
   r = b - Ax;
   r_nrm(1) = norm( r );
%    r_nrm(1) = norm( r )/norm( b );

   if ( r_nrm(1) <= tol ) 
       return;
   else
   end
   iter = 1;
   notconv = 1;
   
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
        w = mats*V(:,i);
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
          
      y = H(1:i,1:i) \ s(1:i);
      x = x + V(:,1:i)*y;
   iter = iter - 1;
   r_nrm = r_nrm(1:iter); % set to actual size
end