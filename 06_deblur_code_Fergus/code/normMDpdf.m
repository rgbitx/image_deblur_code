function y=normMDpdf(x,mu,sig)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

%
% x should be nDims x nPoints 
% mu should be nDims x 1
% sig should be nDims x nDims
%
% y is 1 x nPoints

%if size(x,1)>size(x,2)
%   x=x';
%end

nPoints = size(x,2);
nDims   = size(x,1);

mu_block= mu * ones(1,nPoints); 
i_sig=inv(sig);
d=((2*pi)^-(nDims/2))/sqrt(det(sig));

tt=x-mu_block;
ttt=i_sig*tt;
e=sum(tt.*ttt);

y = d*exp(-0.5*e);




 
