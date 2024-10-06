function [x,J] = mirls(x0,b, Param, mu, max_iter_cgls, max_iter_irls, tol1,tol2)

Wr = ones(size(b));
Wx = ones(size(x0));
u0 = x0;

 Diff = 99999.;

 for j = 1 : max_iter_irls

     u = mcglsw(u0, b ,Param, Wr, Wx, mu, max_iter_cgls, tol1, 0);
     x = Wx.*u;
     Wx = sqrt(abs(x) + 0.00001);  
     e = Wr.*(operator(x,Param,1)-b);
     
     u0=x;
     
     J(j) = sum(abs(e(:)).^2)+mu*sum(abs(x(:)));

  if j>1; 
    Diff = abs(J(j) - J(j-1))/((J(j)+J(j-1))/2);
  end

  if Diff<tol2
   break
  end

 
 end

%  fprintf(' IRLS ended after  %4.0f iterations of %4.0f  \n',  j,max_iter_irls)

return
