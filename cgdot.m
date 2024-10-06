function [out] = cgdot(in1,in2);

 temp  =   in1.*conj(in2);
 out = sum(temp(:));

return;


