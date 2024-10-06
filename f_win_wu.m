function [ h_win ] = f_win_wu(low1,low2,high1,high2,N,dt)
df=1/dt/N;
h_win=zeros(1,N);
ilow1=floor(low1/df);
ilow2=floor(low2/df);
ihigh1=floor(high1/df);
ihigh2=floor(high2/df);
                  
k1=ilow1:ilow2;
k2=ihigh1:ihigh2;

if ilow1~=ilow2
    h_win(k1+1)=cos(pi/2*(ilow2-k1)/(ilow2-ilow1));
    h_win(ilow2+1:ihigh1-1)=1;
    h_win(k2)=cos(pi/2*(k2-ihigh1)/(ihigh2-ihigh1));
end
if ilow1==ilow2
    h_win(ilow1+1:ihigh1-1)=1;
    h_win(k2)=cos(pi/2*(k2-ihigh1)/(ihigh2-ihigh1));
end

if mod(N,2)==0; 
    h_win(N/2+2:N)=fliplr(h_win(2:N/2));
else 
    h_win((N+1)/2+1:N)=fliplr(h_win(2:(N+1)/2));
end  
end

