function B_matrix=Blocking_matrix_f_win_3(C1,C2,C3)
[rowS,colS]=size(C1);
B_matrix=zeros(3*rowS,3*colS);
for i=1:3
    if i==1
        C=C1;
    elseif i==2
        C=C2;
    else
        C=C3;
    end
    B_matrix((i-1)*rowS+1:i*rowS,colS*(i-1)+1:colS*i)=C;
    C=[];
end
end
