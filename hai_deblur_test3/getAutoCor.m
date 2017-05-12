function A=getAutoCor(x,k_sz1,k_sz2)

%this function provides an extra efficent way to compute the
%autocorr of all [k_sz1xk_sz2] windows in x.
%It is equivalent to computing the following
%xM=im2col(x,[k_sz1,k_sz2],'sliding'); %extract all [k_sz1xk_sz2] windows
%A=A+xM*xM';

[M1,M2]=size(x);
sM1=M1-k_sz1+1;
sM2=M2-k_sz2+1;

    
k_sz=k_sz1*k_sz2;
A=zeros(k_sz,k_sz);



for d2=0:k_sz2-1
  for d1=0:k_sz1-1
   
    xx=x(1:end-d1,1:end-d2).*x(d1+1:end,d2+1:end);
    cs=cumsum(cumsum(xx,1),2);

    for j2=0:k_sz2-1
      for j1=0:k_sz1-1
         i1=j1+d1; i2=j2+d2;
         i=(i2)*k_sz1+i1+1;
         j=(j2)*k_sz1+j1+1;
         if ((i>k_sz)|(i1>=k_sz1)|(i2>=k_sz2)|(i1<0)|(i2<0))
            continue
         end

         ts=cs(j1+sM1,j2+sM2);
         if (j1>0)
           ts=ts-cs(j1,j2+sM2);
         end
         if (j2>0)
           ts=ts-cs(j1+sM1,j2);
         end
         if ((j1>0)&(j2>0))
           ts=ts+cs(j1,j2);
         end
         A(j,i)=ts;
         A(i,j)=ts;
      end
    end

  end
end




for d2=0:k_sz2-1
  for d1=1:k_sz1-1

    xx=x(d1+1:end,1:end-d2).*x(1:end-d1,d2+1:end);
    cs=cumsum(cumsum(xx,1),2);
    
    for j2=0:k_sz2-1
      for j1=d1:k_sz1-1
         i1=j1-d1; i2=j2+d2;
         i=(i2)*k_sz1+i1+1;
         j=(j2)*k_sz1+j1+1;
         if ((i>k_sz)|(i1>=k_sz1)|(i2>=k_sz2)|(i1<0)|(i2<0))
            continue
         end
       
         ts=cs(j1-d1+sM1,j2+sM2);
         if (j1>d1)
           ts=ts-cs(j1-d1,j2+sM2);
         end
         if (j2>0)
           ts=ts-cs(j1-d1+sM1,j2);
         end
         if ((j1>d1)&(j2>0))
           ts=ts+cs(j1-d1,j2);
         end
         A(j,i)=ts;
         A(i,j)=ts;
      end
    end
   
  end
end
