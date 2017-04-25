function  pos_arr   =  exemplar_matching(im, par)
%S         =   41;
f         =   par.win;
f2        =   f^2;
s         =   par.step;

N         =   size(im,1)-f+1;
M         =   size(im,2)-f+1;
r         =   [1:s:N];
r         =   [r r(end)+1:N];
c         =   [1:s:M];
c         =   [c c(end)+1:M];
L         =   N*M;
X         =   zeros(f*f, L, 'single');

k    =  0;
for i  = 1:f
    for j  = 1:f
        k    =  k+1;
        blk  =  im(i:end-f+i,j:end-f+j);
        X(k,:) =  blk(:)';
    end
end

% Index image
I     =   (1:L);
I     =   reshape(I, N, M);
N1    =   length(r);
M1    =   length(c);
pos_arr   =  zeros(par.nblk, N1*M1 );
X         =  X';

for  i  =  1 : N1
    for  j  =  1 : M1
        
        row     =   r(i);
        col     =   c(j);
        off     =  (col-1)*N + row;
        off1    =  (j-1)*N1 + i;
                
        template=X(off,:);
        temp=repmat(template,L,1);
        dif=(temp-X)';
        dis=mean(dif.*dif);
        [val,ind]   =  sort(dis);        
        pos_arr(:,off1)  =  ind(1:par.nblk) ;        
    end
end