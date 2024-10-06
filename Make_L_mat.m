function [m]=Make_L_mat(p)

m=p;
Dv = spdiags([-ones(m,1) ones(m,1)],[0 1],m-1,m);
Dv(m,:) = 0;
m=Dv;
end