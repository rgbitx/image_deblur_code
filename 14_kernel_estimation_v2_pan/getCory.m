function b=getCory(x,y,k_sz1,k_sz2)

[n1,n2]=size(x);
sn1=n1-k_sz1+1;
sn2=n2-k_sz2+1;

k_sz=k_sz1*k_sz2;
b=zeros(k_sz,1);
d=0;
for d2=0:k_sz2-1
  for d1=0:k_sz1-1
    d=d+1;
    b(d)=sum(sum(x(1+d1:end-k_sz1+1+d1,1+d2:end-k_sz2+1+d2).*y));
  end
end
