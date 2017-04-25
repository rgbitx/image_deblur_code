function  X  =  Im2Patch( im, par )
f       =   par.win;
N       =   size(im,1)-f+1;
M       =   size(im,2)-f+1;
L       =   N*M;
X       =   zeros(f*f, L, 'single');
k       =   0;
for i  = 1:f
    for j  = 1:f
        k    =  k+1;
        blk  =  im(i:end-f+i,j:end-f+j);
        X(k,:) =  blk(:)';
    end
end