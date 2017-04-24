function [ ] = ultimateSubplot ( NX, NY, x, y, marg )

% Author: Markus Weber
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology
 
if isempty(y),
    
    n = x;
    
    y = ceil(n / NX);
    
    x = n - (y - 1) * NX;
end 

widX = 1 / (NX + (NX + 1) * marg);
widY = 1 / (NY + (NY + 1) * marg);

subplot('position', [(x - 1)*widX + x*widX*marg, (NY - y) * widY + (NY - y + 1) * widY * marg, widX, widY]);
     
     

