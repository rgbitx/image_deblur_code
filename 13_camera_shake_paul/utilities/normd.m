function B = normd(A,p,dim)
% normd(A,p,dim) computes the norm of the subvectors in array A
% along the dimension dim.
% EXAMPLE: if B = normd(A,2,3) then B(i,j,:) = norm(A(i,j,:),2).

B = sum(abs(A).^p,dim).^(1/p);