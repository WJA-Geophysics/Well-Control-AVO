function V=decorrelation(X,iSample_num)
% X=[Rvp Rvs Rrho];
avr=mean(X,1);
aX=[ones(iSample_num,1)*avr(1),ones(iSample_num,1)*avr(2),ones(iSample_num,1)*avr(2)];
Cr=(X-aX)'*(X-aX)./(iSample_num-1);
Cr=X'*X./(iSample_num-1); 

[U,S,V]=svd(Cr); 
v=S^(-1/2)*U^-1;
V=kron(v,eye(iSample_num,iSample_num));
V=V^-1;

