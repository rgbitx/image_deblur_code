function [x]=deconv_ansio_L1(I,k,we,max_it)
% TV deblurring according to 
%Levin et al's SIGGRAPH07 paper.

if (~exist('max_it','var'))
   max_it  =200;
end

[N1,N2]=size(I);

hfs1_x1=floor((size(k,2)-1)/2);
hfs1_x2=ceil((size(k,2)-1)/2);
hfs1_y1=floor((size(k,1)-1)/2);
hfs1_y2=ceil((size(k,1)-1)/2);

hfs_x1=hfs1_x1;
hfs_x2=hfs1_x2;
hfs_y1=hfs1_y1;
hfs_y2=hfs1_y2;


N2=N2+hfs_x1+hfs_x2;
N1=N1+hfs_y1+hfs_y2;


tI=I;
I=zeros(N1,N2);
I(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=tI; 
x=I;


dxf=[1 -1];
dyf=[1;-1];


weight_x=ones(N1,N2-1);
weight_y=ones(N1-1,N2);


[x]=deconvL1_irls(x(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2),k,we,max_it,weight_x,weight_y);


w0=0.1;
exp_a=1;
thr_e=0.01; 

for t=1:2

  dy=conv2(x,fliplr(flipud(dyf)),'valid');
  dx=conv2(x,fliplr(flipud(dxf)),'valid');


  weight_x=w0*max(abs(dx),thr_e).^(exp_a-2); 
  weight_y=w0*max(abs(dy),thr_e).^(exp_a-2);
  
  [x]=deconvL1_irls(I(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2),k,we,max_it,weight_x,weight_y);

end


x=x(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2);

return


