function cI=fftconvf(I,k,K,method)

[N1,N2]=size(I);
[k1,k2]=size(k);
hk1=(k1-1)/2;
hk2=(k2-1)/2;
[bk1,bk2]=size(K);
hdiff1d=ceil((bk1-N1)/2);
hdiff1u=floor((bk1-N1)/2);
hdiff2d=ceil((bk2-N2)/2);
hdiff2u=floor((bk2-N2)/2);


I=zero_pad2(I,hdiff1d,hdiff1u,hdiff2d,hdiff2u);

fI=fft2(ifftshift(I));


cI=fftshift(ifft2(fI.*K));

if exist('method','var')

if strcmp(method,'same')
  cI=cI(hdiff1d+1:end-hdiff1u,hdiff2d+1:end-hdiff2u);    
end
if strcmp(method,'valid')
  cI=cI(hdiff1d+1:end-hdiff1u,hdiff2d+1:end-hdiff2u);    
  cI=cI(hk1+1:end-hk1,hk2+1:end-hk2);     
end

end
