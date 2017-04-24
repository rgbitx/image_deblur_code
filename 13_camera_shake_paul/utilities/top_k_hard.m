function x = top_k_hard(y,k)

sy = sign(y);
ay = sort(abs(y(y ~= 0)),'descend');
x = sy.*y.*(abs(y) >= ay(k));

end