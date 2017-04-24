% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

lambda = [0.1 2 8];
pi = [0.3 0.3 0.3];

pi = pi / sum(pi);

x = [0:0.01:10];
y = zeros(size(x));

figure(100); clf;

for a=1:length(lambda)
  
  y = y + pi(a) * lambda(a) * exp(- lambda(a) * x);
  
  yy(a,:) =  lambda(a) * exp(- lambda(a) * x);

  plot(x,yy(a,:),'k--');
  hold on;
end

plot(x,y,'k');





