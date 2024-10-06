function G = operaterG(wt, M)
% wt: seismic wavelet
% M: sampling of seismic data
[nt,nx] = size(wt);
if nt<nx
    wt = wt';
end
nw = length(wt);
G = convmtx(wt,M);
n1 = floor(nw/2);
z=n1+1:M+nw-1-n1;
if length(z)==M
    G = G(n1+1:M+nw-1-n1,:);
else
    G = G(n1:M+nw-1-n1,:);
end

end