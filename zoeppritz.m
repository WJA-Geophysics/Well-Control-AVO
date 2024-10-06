function coef=zoeppritz(rho1,a1,b1,rho2,a2,b2,incwav,irfwav,ipol,anginc)		
rpd = pi/180.0;
%
% ======================= Solid-Solid interface ========================
if (a1&a2&b1&b2)>0
      if incwav==1
         i1 = anginc * rpd;
         p = sin(i1)/a1;
         ci1 = cos(i1);
         cj1 = sqrt(1. - (p*b1).^2);
      elseif incwav==2
         j1 = anginc * rpd;
         p = sin(j1)/b1;
         cj1 = cos(j1);
         ci1 = sqrt(1. - (p*a1).^2);
      end
%
      ci2 = sqrt(1. - (p*a2).^2);
      cj2 = sqrt(1. - (p*b2).^2);
      ca1 = ci1/a1;
      ca2 = ci2/a2;
      cb1 = cj1/b1;
      cb2 = cj2/b2;
      rb1 = rho1 * (1. - 2.*(b1*p).^2);
      rb2 = rho2 * (1. - 2.*(b2*p).^2);
      a = rb2 - rb1;
      b = rb2 + 2.*rho1*(b1*p).^2;
      c = rb1 + 2.*rho2*(b2*p).^2;
      d = 2.*(rho2*b2^2 - rho1*b1.^2);
      e = b.*ca1 + c.*ca2;
      f = b.*cb1 + c.*cb2;
      g = a - d.*ca1.*cb2;
      h = a - d.*ca2.*cb1;
      dd = e.*f + g.*h.*p.*p;
%
      if incwav==1
         if irfwav==1
            coef = ((b.*ca1 - c.*ca2).*f - (a + d.*ca1.*cb2).*h.*p.*p)./dd;
         elseif irfwav==2
            coef = -(2*ca1.*(a.*b + c.*d.*ca2.*cb2).*p.*a1)./(b1.*dd);
         elseif irfwav==3
            coef = (2*rho1.*ca1.*f.*a1)./(a2.*dd);
         elseif irfwav==4
            coef = (2*rho1.*ca1.*h.*p.*a1)./(b2.*dd);
         end
      elseif incwav==2
         if irfwav==1
            coef = -(2*cb1.*(a.*b + c.*d.*ca2.*cb2).*p.*b1)./(a1.*dd);
         elseif irfwav==2
            coef = -((b.*cb1 - c.*cb2).*e - (a + d.*ca2.*cb1).*g.*p.*p)./dd;
         elseif irfwav==3
            coef = -(2*rho1.*cb1.*g.*p.*b1)./(a2.*dd);
         elseif irfwav==4
            coef =  (2*rho1.*cb1.*e.*b1)./(b2.*dd);
         end
      end
%
      if ipol==1
         ampl = sqrt(real(coef).^2 + imag(coef).^2);
         if (coef==0.)
            phas = 0.;
         else
            phas = atan2(imag(coef), real(coef));
         end
         coef = ampl + i*phas;
      end
%
      return
	end
%
%
% =============== The case of a liquid-solid interface ====================
%
if (a1>0)&(a2>0)&(b1==0)&(b2>0)
   if incwav==1
      i1 = anginc * rpd;
      p = sin(i1)/a1;
      ci1 = cos(i1);
      ci2 = sqrt(1. - (p*a2).^2);
      cj2 = sqrt(1. - (p*b2).^2);
      c2j2 = 1. - 2.*(b2*p).^2;
      t1 = rho1*a1*ci2;
      t2 = rho2*a2*ci1.*c2j2.^2;
      t3 = 4.*rho2*b2^3*p.*p.*ci1.*ci2.*cj2;
      dd = t1 + t2 + t3;
%
      if irfwav==1
         coef = (-t1 + t2 + t3)./dd;
      elseif irfwav==2
         error('irfwav cannot equal 2, execution stopped')
      elseif irfwav==3
         coef = (2*rho1*a1*ci1.*c2j2)./dd;
      elseif irfwav==4
         coef = -(4*rho1*a1*b2*p.*ci1.*ci2)./dd;
      end
%
      if ipol==1
         ampl = sqrt(real(coef).^2 + imag(coef).^2);
         if coef==0
            phas = 0.;
         else
            phas = atan2(imag(coef), real(coef));
         end
         coef = ampl+i*phas;
      end
%  
      return
    elseif incwav==2
    error('You must have an incident P wave');
    end 
end
%
% =============== The case of a solid-liquid interface ===================
%  
if (a1>0)&(a2>0)&(b1>0)&(b2==0)
      if incwav==1
         i1 = anginc * rpd;
         p = sin(i1)/a1;
         ci1 = cos(i1);
         cj1 = sqrt(1. - (p*b1).^2);
      elseif incwav==2
         j1 = anginc * rpd;
         p = sin(j1)/b1;
         cj1 = cos(j1);
         ci1 = sqrt(1. - (p*a1).^2);
      end
%
      ci2 = sqrt(1. - (p*a2).^2);
      c1j1 = 1. - 2.*(b1*p).^2;
      t1 = rho2*a2*ci1;
      t2 = rho1*a1*ci2.*c1j1.^2;
      t3 = 4*rho1*b1^3*p.*p.*ci2.*ci1.*cj1;
      dd = t1 + t2 + t3;
%
      if incwav==1
         if irfwav==1
            coef = (t1 - t2 + t3)./dd;
         elseif irfwav==2
            coef = (4*rho1*a1*b1*p.*c1j1.*ci1.*ci2)./dd;
         elseif irfwav==3
            coef = (2*rho1*a1*c1j1.*ci1)./dd;
         elseif irfwav==4
            error('irfwav cannot equal 4, execution stopped')
         end
      elseif incwav==2
         if irfwav==1
            coef = (4*rho1*b1^2*p.*c1j1.*ci2.*cj1)./dd;
         elseif irfwav==2
            coef = (t1 + t2 - t3)./dd;
         elseif irfwav==3
            coef = -(4*rho1*b1^2*p.*ci1.*cj1)./dd;
         elseif irfwav==4
            error('irfwav cannot equal 4, execution stopped')
         end
      end
%
      if ipol==1
         ampl = sqrt(real(coef).^2 + imag(coef).^2);
         if coef==0.
            phas = 0.;
         else
            phas = atan2(imag(coef), real(coef));
         end
         coef = ampl+i*phas;
      end
%
      return
	end
%
% ================= The case of a liquid-liquid interface =================

if (a1>0)&(a2>0)&(b1==0)&(b2==0)
  if incwav==1
      i1 = anginc * rpd;
      p = sin(i1)/a1;
      ci1 = cos(i1);
      ci2 = sqrt(1. - (p*a2).^2);
      ra1 = rho1*a1*ci2;
      ra2 = rho2*a2*ci1;
      dd = ra1 + ra2;
      if irfwav==1
         coef = (ra2 - ra1)./dd;
      elseif irfwav==2
         error('irfwav cannot equal 2, execution stopped')
      elseif irfwav==3
         coef = (2*rho1*a1*ci1)./dd;
      elseif irfwav==4
         error('irfwav cannot equal 4, execution stopped')
      end
%
      if ipol==1 
         ampl = sqrt(real(coef).^2 + imag(coef).^2);
         if (coef==0.)
            phas = 0.;
         else
            phas = atan2(imag(coef), real(coef));
         end
         coef = ampl+i*phas;
      end
%
      return
   else
   error('You must have an incident P wave')
   end
end
      

%
% =============== The case of SH reflection and transmission ==============
%  
if (a1==00)&(a2==00)&(b1>0)&(b2>0)
   if incwav==2
      j1=anginc*rpd;
      p=sin(j1)/b1;
      cj1=cos(j1);
      cj2=sqrt(1-(p*b2).^2);
      rs1=rho1*b1*cj1;
      rs2=rho2*b2*cj2;
      dd=rs1+rs2;
      if irfwav==1 
         error('irfwav cannot equal 1, execution stopped')
      elseif irfwav==2
         coef=(rs1-rs2)./dd;
      elseif irfwav==3
         error('irfwav cannot equal 3, execution stopped')
      elseif irfwav==4
         coef=(2*rs1)./dd;
      end

      if ipol==1
         ampl=sqrt(real(coef).^2+imag(coef).^2);
         if coef==0
            phas=0;
         else
            phas=atan2(imag(coef),real(coef));
         end
         coef=ampl+i*phas;
      end
   elseif incwav==1
   error('You must have an incident S wave');
   end
end
