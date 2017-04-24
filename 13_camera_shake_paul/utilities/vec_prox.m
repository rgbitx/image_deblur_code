function X = vec_prox(prox,Y,tau)
   nY = normd(Y,2,2);
   nX = prox(nY,tau);
   X = bsxfun(@times,Y,nX./nY);
   X(isnan(X)) = 0;
end