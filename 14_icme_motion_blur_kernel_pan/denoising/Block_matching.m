function  pos_arr   =  Block_matching(im, par)
S         =   21;
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
                
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, N );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, M );
         
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        B       =   X(idx, :);        
        v       =   X(off, :);
        
        
        dis     =   (B(:,1) - v(1)).^2;
        for k = 2:f2
            dis   =  dis + (B(:,k) - v(k)).^2;
        end
        dis   =  dis./f2;
        [val,ind]   =  sort(dis);        
        pos_arr(:,off1)  =  idx( ind(1:par.nblk) );        
    end
end