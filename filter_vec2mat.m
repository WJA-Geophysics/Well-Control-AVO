function  Fil_matrix=filter_vec2mat(fil_vec)

length=numel(fil_vec);
Fil_matrix=zeros(length,length);
for i=1:length
    Fil_matrix(i,i)=fil_vec(i);
end
end