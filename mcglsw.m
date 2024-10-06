function [x] = mcglsw(x0, b ,Param, Wr, Wx, mu, max_iter, tol, prnt)

 b = Wr.*b;

   x = x0;
   r   = b - Wr.*operator(Wx.*x,Param,1);
   s   = Wx.*operator(Wr.*r,Param,-1) - mu*x;
   p = s;

   gamma  = cgdot(s,s);
   norms0 = sqrt(gamma);           % norm of the gradient is used to stop 
   k      = 0;
   flag   = 0;

if prnt 
    
    fprintf( ' ============================================== \n');
    fprintf( ' ================= CGLS ======================= \n');
    fprintf( ' ============================================== \n');
    
    head = '     k           |grad|       |grad|/|grad_0|        '; 
    form = '   %3.0f       %12.5g           %8.3g              \n';
    disp('  ');   disp(head);
    fprintf(form, k, norms0,1)

end

 
 while (k <max_iter) && (flag==0);
  
    q = Wr.*operator(Wx.*p,Param,1);
    delta = cgdot(q,q) + mu*cgdot(p,p);
    if delta == 0, delta = 1.e-10; end
    alpha = gamma/delta; 
    x = x + alpha*p;
    r = r - alpha*q;
    s = Wx.*operator(Wr.*r,Param,-1) - mu*x;

    gamma1 = cgdot(s,s);
    norms = sqrt(gamma1);
    beta = gamma1/gamma;
    gamma = gamma1;
    
    p = s + beta*p;
    
    flag = (norms<=norms0 * tol);
    nres = norms / norms0;

    k = k+1;

 if prnt, fprintf(form, k, norms, nres); end
 
 end; 

 
 % Diagnostics
 
 if  k == max_iter; flag = 2; end

if prnt,
     

               fprintf( ' ============================================== \n');
   if flag==1; fprintf( ' ====== CGLS converged before max_iter ======== \n'); end
   if flag==2; fprintf( ' ====== CGLS reached max_iter ================= \n'); end
               fprintf( ' ============================================== \n \n \n');

end

return