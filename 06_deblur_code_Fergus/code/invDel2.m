function [invK]=invDel2(isize);

% Author: Yair Weiss
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology


 K=zeros(isize);
 K(isize/2,isize/2)=-4;
 K(isize/2+1,isize/2)=1;
 K(isize/2,isize/2+1)=1;
 K(isize/2-1,isize/2)=1;
 K(isize/2,isize/2-1)=1;
 
 Khat=fft2(K);
 I=find(Khat==0);
 Khat(I)=1;
 invKhat=1./Khat;
 invKhat(I)=0;
 invK=ifft2(invKhat);
 invK=-real(invK);
 invK=conv2(invK,[1 0 0;0 0 0;0 0 0],'same');% shift by one
 
