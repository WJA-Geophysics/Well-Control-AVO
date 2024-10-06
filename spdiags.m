function [res1,res2] = spdiags(arg1,arg2,arg3,arg4)

if nargin <= 2
    % Extract diagonals
    A = arg1;
    if nargin == 1
        % Find all nonzero diagonals
        [i,j] = find(A);
        % Compute d = unique(d) without extra function call
        d = sort(j-i);
        d = d(diff([-inf; d(:)])~=0);
        d = d(:);
    else
        % Diagonals are specified
        d = arg2(:);
    end
    [m,n] = size(A);
    p = length(d);
    B = zeros(min(m,n),p,class(A));
    for k = 1:p
        if m >= n
            i = max(1,1+d(k)):min(n,m+d(k));
        else
            i = max(1,1-d(k)):min(m,n-d(k));
        end
        B(i,k) = diagk(A,d(k));
    end
    res1 = B;
    res2 = d;
end

if nargin >= 3
    B = arg1;
    d = arg2(:);
    p = length(d);
    if nargin == 3 % Replace specified diagonals
        A = arg3;
    else           % Create new matrix with specified diagonals
        A = sparse(arg3, arg4);
    end
    [m,n] = size(A);
    
    % Check size of B. Should be min(m,n)-by-p for matrix,
    % min(m,n) for column vector, p for row vector, or a scalar.
    % For backwards compatibility, only error if the code would
    % previously have errored out in the indexing expression.
    [mB, nB] = size(B);
    maxIndexRows = max(max(1,1-d), min(m,n-d)) + (m>=n)*d;
    maxIndexRows(max(1,1-d) > min(m,n-d)) = 0;
    if (mB ~= 1 && any(maxIndexRows > mB)) || (nB ~= 1 && nB < p)
        if nargin == 3
            error(message('MATLAB:spdiags:InvalidSizeBThreeInput'));
        else
            error(message('MATLAB:spdiags:InvalidSizeBFourInput'));
        end
    end

    % Compute indices and values of sparse matrix with given diagonals
    if nnz(A) == 0 && (~isobject(B) && ~isobject(d) && ~isobject(m) && ~isobject(n)) && ...
            (mB == 1 || mB == min(m, n)) && (nB == 1 || nB == p)
        B = full(B);
        
        if m < n
            % Compute transpose of A, then transpose before returning.
            d = -d;
        end
        
        % Sort d in descending order and reorder B accordingly:
        if issorted(d, 'descend')
        elseif issorted(d, 'ascend')
            d = flip(d);
            B = flip(B, 2);
        else
            [d, ind] = sort(d, 'descend');
            if ~iscolumn(B)
                B = B(:, ind);
            end
        end

        % Merge duplicate diagonals
        hasDuplicates = nnz(diff(d)) < length(d)-1;
        if hasDuplicates
            [B, d] = mergeDuplicateDiagonals(B, d);
        end
        
        if m >= n
            res1 = matlab.internal.sparse.makeSparsePresorted(double(B), double(d), m, n);
        else
            % Compute transpose of A, then transpose before returning.
            res1 = matlab.internal.sparse.makeSparsePresorted(double(B), double(d), n, m);
            res1 = res1.';
        end
    else
        % Insert diagonals into existing matrix A. This algorithm
        % has fewer requirements on B than the one above.
        res1 = makeSparseGeneral(B, d, A);
    end
    
    if islogical(A) || islogical(B)
        res1 = (res1~=0);
    end
end
end


function A = makeSparseGeneral(B, d, A)
% Construct sparse matrix by inserting
% diagonals into existing matrix A.

[m, n] = size(A);
p = length(d);

% Precompute number of nonzeros (parts of each diagonal that overlap with
% the matrix) and allocate inputs I, J and V for sparse
nz = sum(max(0, min(m, n-d) - max(1, 1-d) + 1));
I = zeros(nz, 1);
J = zeros(nz, 1);
V = zeros(nz, 1);

% Fill in the diagonals
offset = 1;
if isscalar(B)
    for k = 1:p
        % Append new d(k)-th diagonal to compact form
        for i=max(1, 1-d(k)):min(m, n-d(k))
            I(offset) = i;
            J(offset) = i + d(k);
            V(offset) = B;
            offset = offset + 1;
        end
    end
elseif isrow(B)
    for k = 1:p
        % Append new d(k)-th diagonal to compact form
        for i=max(1, 1-d(k)):min(m, n-d(k))
            I(offset) = i;
            J(offset) = i + d(k);
            V(offset) = B(k);
            offset = offset + 1;
        end
    end
elseif iscolumn(B)
    for k = 1:p
        % Append new d(k)-th diagonal to compact form
        for i=max(1, 1-d(k)):min(m, n-d(k))
            I(offset) = i;
            J(offset) = i + d(k);
            if m >= n
                V(offset) = B(i + d(k));
            else
                V(offset) = B(i);
            end
            offset = offset + 1;
        end
    end
else
    for k = 1:p
        % Append new d(k)-th diagonal to compact form
        for i=max(1, 1-d(k)):min(m, n-d(k))
            I(offset) = i;
            J(offset) = i + d(k);
            if m >= n
                V(offset) = B(i + d(k), k);
            else
                V(offset) = B(i, k);
            end
            offset = offset + 1;
        end
    end
end

if nnz(A) > 0
    % Process A in compact form
    [Iold,Jold,Vold] = find(A);
    
    % Delete current d(k)-th diagonal, k=1,...,p
    i = any((Jold(:) - Iold(:)) == d', 2);
    Iold(i) = [];
    Jold(i) = [];
    Vold(i) = [];
    
    % Combine new diagonals and non-diagonal entries of original matrix
    I = [I(:); Iold(:)];
    J = [J(:); Jold(:)];
    V = [V(:); Vold(:)];
end

A = sparse(I, J, V, m, n);
end


function D = diagk(X,k)
% DIAGK  K-th matrix diagonal.
% DIAGK(X,k) is the k-th diagonal of X, even if X is a vector.
D = matlab.internal.math.diagExtract(X,k);
end


function [BinMerged, dMerged] = mergeDuplicateDiagonals(Bin, d)
% Combine columns addressed to the same diagonal
if ~iscolumn(Bin)
    [dMerged, ~, input2output] = unique(d, 'stable');
    M = sparse(1:length(d), input2output, 1);
    BinMerged = Bin * M;
else
    [dMerged, ~, input2output] = unique(d, 'stable');
    repetitions = accumarray(input2output(:), 1);
    BinMerged = Bin .* repetitions';
end
end