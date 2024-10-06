function out = operator(in,Param,flag);
% wrapper to use cgls.m when the linear operator is a matrix 

A = Param;

if flag==1; 
    out = A*in; 
else;
    out = A'*in;
end 
