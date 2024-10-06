function [x]=Xi_matrix(model,p)
x=zeros([p,1]);
for i=1:p
    x(i,1)=0.5*log(model(i,1)/model(1,1));%);
end
end