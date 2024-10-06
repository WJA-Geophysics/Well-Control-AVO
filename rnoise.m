function noise=rnoise(v,s_to_n,irange,flag)

 if nargin<4
   flag=1;
 end
 if nargin<3
   irange=1:length(v);
 end
 if nargin<2
   s_to_n=2;
 end
% 
 done=0;
 if flag==1
  done=1;
 end
 if flag==0
   done=1; 
 end
 if done==0
   error(' invalid flag');
 end
%
 c=clock;
 n=fix(10.*c(6));
 if flag == 1
    randn('seed',n);
 end
 if flag == 0
	rand('seed',n);
  end  
 noise=rand(size(v));
% adjust to zero mean
 noise=noise-mean(noise);

% measure signal and noise powers
 [m,n]=size(v);
 if m==1
   ps= sqrt(v(irange)*v(irange)');
   pn= sqrt(noise(irange)*noise(irange)');
   scalar1= ps/(pn*s_to_n);
 end
 if n==1 
   ps= sqrt(v(irange)'*v(irange));
   pn= sqrt(noise(irange)'*noise(irange));
   scalar1= ps/(pn*s_to_n);
 end
% adjust noise power
 noise=noise*scalar1;


 
