function [x]=deconvL1_irls(I,k,we,max_it,weight_x,weight_y)

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

if (max(size(k)<=5))
  Ax=conv2(conv2(x,fliplr(flipud(k)),'same').*mask,  k,'same');
else
  Ax=fftconvf(fftconvf(x,fliplr(flipud(k)),conj(K),'same').*mask,k,K,'same');
end


Ax=Ax+we*conv2(weight_x.*conv2(x,fliplr(flipud(dxf)),'valid'),dxf);
Ax=Ax+we*conv2(weight_y.*conv2(x,fliplr(flipud(dyf)),'valid'),dyf);


r = b - Ax;

for iter = 1:max_it  
     rho = (r(:)'*r(:));
     if (rho<0.1^8)
       %iter
       %'convarged'
       break
     end
     
     if ( iter > 1 ),                      
        beta = rho / rho_1;
        p = r + beta*p;
     else
        p = r;
     end
     if (max(size(k)<5))
       Ap=conv2(conv2(p,fliplr(flipud(k)),'same').*mask,  k,'same');
     else  
       Ap=fftconvf(fftconvf(p,fliplr(flipud(k)),conj(K),'same').*mask,k,K,'same');
     end

     Ap=Ap+we*conv2(weight_x.*conv2(p,fliplr(flipud(dxf)),'valid'),dxf);
     Ap=Ap+we*conv2(weight_y.*conv2(p,fliplr(flipud(dyf)),'valid'),dyf);


     q = Ap;
     alpha = rho / (p(:)'*q(:) );
     x = x + alpha * p;                    

     r = r - alpha*q;                     

     rho_1 = rho;
end
end
%%
function N=goodfactor(N)
    
    f=factor(N);
    
    while(max(f)>7)
      N=N+1;
      f=factor(N);
    end
end



