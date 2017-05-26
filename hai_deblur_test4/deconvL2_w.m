function [x]=deconvL2_w(I,k,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy)

if (~exist('max_it','var'))
   max_it=200;
end

[N1,N2]=size(I);



[fs_y,fs_x]=size(k);
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
N=N2*N1;
mask=zeros(N1,N2);
mask(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=1;

if (~exist('weight_x','var'))
  weight_x=ones(N1,N2-1);
  weight_y=ones(N1-1,N2);
  weight_xx=zeros(N1,N2-2);
  weight_yy=zeros(N1-2,N2);
  weight_xy=zeros(n-1,m-1);
end


tI=I;
I=zeros(N1,N2);
I(hfs_y1+1:N1-hfs_y2,hfs_x1+1:N2-hfs_x2)=tI; 
x=tI([ones(1,hfs_y1),1:end,end*ones(1,hfs_y2)],[ones(1,hfs_x1),1:end,end*ones(1,hfs_x2)]);



b=conv2(x.*mask,k,'same');

%pad k with zeros up to a nearby integer with small prime factors, for fast fft
N1p=goodfactor(N1+hfs1_y1+hfs1_y2);
N2p=goodfactor(N2+hfs1_x1+hfs1_x2);
K=zero_pad2(k,ceil((N1p-fs_y)/2),floor((N1p-fs_y)/2),ceil((N2p-fs_x)/2),floor((N2p-fs_x)/2));
K=fft2(ifftshift(K));



dxf=[1 -1];
dyf=[1;-1];
dyyf=[-1; 2; -1];
dxxf=[-1, 2, -1];
dxyf=[-1 1;1 -1];

if (max(size(k)<=5))
  Ax=conv2(conv2(x,fliplr(flipud(k)),'same').*mask,  k,'same');
else
  Ax=fftconvf(fftconvf(x,flp(k),conj(K),'same').*mask,k,K,'same');
end


Ax=Ax+we*conv2(weight_x.*conv2(x,fliplr(flipud(dxf)),'valid'),dxf);
Ax=Ax+we*conv2(weight_y.*conv2(x,fliplr(flipud(dyf)),'valid'),dyf);
Ax=Ax+we*(conv2(weight_xx.*conv2(x,fliplr(flipud(dxxf)),'valid'),dxxf));
Ax=Ax+we*(conv2(weight_yy.*conv2(x,fliplr(flipud(dyyf)),'valid'),dyyf));
Ax=Ax+we*(conv2(weight_xy.*conv2(x,fliplr(flipud(dxyf)),'valid'),dxyf));


r = b - Ax;

for iter = 1:max_it  
     rho = (r(:)'*r(:));

     if ( iter > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = r + beta*p;
     else
        p = r;
     end
     if (max(size(k)<5))
       Ap=conv2(conv2(p,fliplr(flipud(k)),'same').*mask,  k,'same');
     else  
       Ap=fftconvf(fftconvf(p,flp(k),conj(K),'same').*mask,k,K,'same');
     end

     Ap=Ap+we*conv2(weight_x.*conv2(p,fliplr(flipud(dxf)),'valid'),dxf);
     Ap=Ap+we*conv2(weight_y.*conv2(p,fliplr(flipud(dyf)),'valid'),dyf);
     Ap=Ap+we*(conv2(weight_xx.*conv2(p,fliplr(flipud(dxxf)),'valid'),dxxf));
     Ap=Ap+we*(conv2(weight_yy.*conv2(p,fliplr(flipud(dyyf)),'valid'),dyyf));
     Ap=Ap+we*(conv2(weight_xy.*conv2(p,fliplr(flipud(dxyf)),'valid'),dxyf));


     q = Ap;
     alpha = rho / (p(:)'*q(:) );
     x = x + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual

     rho_1 = rho;
end



