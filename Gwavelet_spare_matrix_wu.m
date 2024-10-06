function Gwavelet_spare = Gwavelet_spare_matrix_wu(Gwavelet,n)
[mw nw]=size(Gwavelet);
Gwavelet_spare=[];

for i=1:1:n
    Gwavelet_spare(mw*(i-1)+1:mw*i,nw*(i-1)+1:nw*i)=Gwavelet;
end

end
