function [x]=deconvSps(I,k,we,max_it)
%This is the original sparse non blind deconvolution code from
%Levin et al's SIGGRAPH07 paper. (In the image space)
%Key differences
%1) For regularization, it uses 2nd order derivatives and not only
%1st order ones 
%2) The prior is simply taken as |z|^0.8, no MOG.
%note: size(k) is expected to be odd in both dimensions 

if (~exist('max_it','var'))
   max_it  =200;
end

[N1,N2]=size(I);

hfs1_x1=floor((size(k,2)-1)/2);
hfs1_x2=ceil((size(k,2)-1)/2);
hfs1_y1=floor((size(k,1)-1)/2);
hfs1_y2=ceil((size(k,1)-1)/2);
shifts1=[-hfs1_x1  hfs1_x2  -hfs1_y1  hfs1_y2];

hfs_x1=hfs1_x1;
hfs_x2=hfs1_x2;
hfs_y1=hfs1_y1;
hfs_y2=hfs1_y2;


N2=N2+hfs_x1+hfs_x2;
N1=N1+hfs_y1+hfs_y2;
N=N1*N2;
mask=zeros(N1,N2);
mask(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=1;


tI=I;
I=zeros(N1,N2);
I(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=tI; 
x=I;


dxf=[1 -1];
dyf=[1;-1];
dyyf=[-1; 2; -1];
dxxf=[-1, 2, -1];
dxyf=[-1 1;1 -1];


weight_x=ones(N1,N2-1);
weight_y=ones(N1-1,N2);
weight_xx=ones(N1,N2-2);
weight_yy=ones(N1-2,N2);
weight_xy=ones(N1-1,N2-1);


[x]=deconvL2_w(x(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2),k,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy);


w0=0.1;
exp_a=0.8;
thr_e=0.01; 

for t=1:2

  dy=conv2(x,flp(dyf),'valid');
  dx=conv2(x,flp(dxf),'valid');
  dyy=conv2(x,flp(dyyf),'valid');
  dxx=conv2(x,flp(dxxf),'valid');
  dxy=conv2(x,flp(dxyf),'valid');


  weight_x=w0*max(abs(dx),thr_e).^(exp_a-2); 
  weight_y=w0*max(abs(dy),thr_e).^(exp_a-2);
  weight_xx=0.25*w0*max(abs(dxx),thr_e).^(exp_a-2); 
  weight_yy=0.25*w0*max(abs(dyy),thr_e).^(exp_a-2);
  weight_xy=0.25*w0*max(abs(dxy),thr_e).^(exp_a-2);
  
  [x]=deconvL2_w(I(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2),k,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy);

end


x=x(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2);

return


