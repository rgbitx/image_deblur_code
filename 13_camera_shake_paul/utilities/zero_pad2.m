function zM=zero_pad2(M,zp1d,zp1u,zp2d,zp2u);


[n,m,k]=size(M);

zM=zeros(n+zp1u+zp1d,m+zp2d+zp2u,k);
zM(zp1d+1:end-zp1u,zp2d+1:end-zp2u,:)=M;
