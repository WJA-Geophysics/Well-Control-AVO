function xsz=angleterm_spare_matrix_wu(g,number)
[mg ng]=size(g);
xsz=[];
ang=[];

for i=1:1:mg
    for j=1:1:ng
        ang=g(i,j)*ones(number,1);
        ang=diag(ang);
        xsz(number*(i-1)+1:number*i,number*(j-1)+1:number*j)=ang;
    end
end
end
       
        