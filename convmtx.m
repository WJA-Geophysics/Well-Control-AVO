function A = convmtx(h,n)
%CONVMTX Convolution matrix.
%   A = CONVMTX(H,N) returns the convolution matrix, A, such that
%   convolution of H and another vector, X, of length N with the same
%   vector orientation as H may be expressed by matrix multiplication:
%   When H and X are column vectors, A*X is the same as CONV(H,X);
%   when H and X are row vectors, X*A is the same as CONV(H,X).
%
%   % Example: 
%   %   Generate a simple convolution matrix.
%
%   h = [1 2 3 2 1];
%   convmtx(h,7)        % Convolution matrix
%
%   See also CONV.

%   Copyright 1988-2016 The MathWorks, Inc.

%#codegen

narginchk(2,2);

validateattributes(h,{'numeric','embedded.fi'},{'vector'},'convmtx','H',1);
validateattributes(n,{'numeric'},{'real','scalar','finite','integer','positive'},'convmtx','N',2);

if coder.target('MATLAB')
    m = length(h) + n - 1;
    A = repmat([h(:); zeros(n,1)], 1, n+1);
    A = reshape(A(1:m*n), m, n);
    
    if isrow(h)
        A = A.';
    end
else
    nh = coder.internal.indexInt(length(h));
    nA = coder.internal.indexInt(n);
    mA = nh + nA - 1;
    if isrow(h)
        A = zeros(nA,mA,'like',h);
        for i = 1:nA
            for j = 1:nh
                A(i,j + i - 1) = h(j);
            end
        end
    else
        A = zeros(mA,nA,'like',h);
        for j = 1:nA
            for i = 1:nh
                A(i + j - 1,j) = h(i);
            end
        end
    end
end
