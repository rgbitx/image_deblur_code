function k=resizeKer(k,ret,k1,k2)
%%
% levin's code
k=imresize(k,ret);
k=max(k,0);
k=fixsize(k,k1,k2);
if max(k(:))>0
    k=k/sum(k(:));
end
end
%% 
function nf=fixsize(f,nk1,nk2)
[k1,k2]=size(f);

while((k1~=nk1)|(k2~=nk2))
    
    if (k1>nk1)
        s=sum(f,2);
        if (s(1)<s(end))
            f=f(2:end,:);
        else
            f=f(1:end-1,:);
        end
    end
    
    if (k1<nk1)
        s=sum(f,2);
        if (s(1)<s(end))
            tf=zeros(k1+1,size(f,2));
            tf(1:k1,:)=f;
            f=tf;
        else
            tf=zeros(k1+1,size(f,2));
            tf(2:k1+1,:)=f;
            f=tf;
        end
    end
    
    if (k2>nk2)
        s=sum(f,1);
        if (s(1)<s(end))
            f=f(:,2:end);
        else
            f=f(:,1:end-1);
        end
    end
    
    if (k2<nk2)
        s=sum(f,1);
        if (s(1)<s(end))
            tf=zeros(size(f,1),k2+1);
            tf(:,1:k2)=f;
            f=tf;
        else
            tf=zeros(size(f,1),k2+1);
            tf(:,2:k2+1)=f;
            f=tf;
        end
    end
    
    
    
    [k1,k2]=size(f);
    
end

nf=f;
end