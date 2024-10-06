function E_matrix=Reverse_extension_matrix(length)
E=eye([length length]);
E_matrix=[fliplr(E);E;fliplr(E)];
end


