function x=randnND(mu,Sigma,N)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

%%  x=randN(mu,Sigma,N)
%% mu, sigma: parameters of the Gaussian
%% N: how many samples should be generated
%% 

   mu = mu(:);
   l = length(mu);
   [V,D] = eig(Sigma); %% Calculate the transformation matrix
   ss = sqrt(diag(D));
   T =  V*diag(ss);
   x = T*randn(l,N) + mu*ones(1,N);
