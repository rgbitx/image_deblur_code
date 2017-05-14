function [psf,opts] = estimate_psf_gradient_mul(B,L,opts,psf_old)
% The objective function: 
% psf^* = argmin ||Lx*k - Bx||^2 +||Ly*k - By||^2+ weight |K|^2
psf_size = size(psf_old);
weight = opts.gamma;
% These values can be pre-computed at the beginning of each level
% FBx,FBy
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];
Lx = conv2(L,dx,'valid');
Ly = conv2(L,dy,'valid');
%[Lx,Ly,opts.grad_threshold]= threshold_pxpy_v1(L,max(psf_size(1:2)),threshold);
FLx = fft2(Lx);
FLy = fft2(Ly);
psf = zeros(psf_size);
Bx = conv2(B(:,:),dx,'valid');
By = conv2(B(:,:),dy,'valid');
FBx = fft2(Bx);
FBy = fft2(By);    
psf(:,:) = estimate_psf(FBx,FBy,FLx,FLy,weight,psf_size(1:2));

end%function 

function psf = estimate_psf(FBx,FBy,FLx,FLy,weight,psf_size)
% compute b = sum_i w_i latent_i * blurred_i
b_f = conj(FLx).*FBx+conj(FLy).* FBy;
b = real(otf2psf(b_f, psf_size));

p.m = conj(FLx).*FLx+conj(FLy).*FLy;
p.img_size = size(FBx);
p.psf_size = psf_size;
p.lambda = weight;

psf = ones(psf_size) / prod(psf_size);
psf = conjgrad(psf,b,20,1e-5, @compute_Ax, p);

% normalized kernel
%psf(psf<0) = 0;
psf(psf<max(psf(:))*0.05) = 0;
psf = psf / sum(psf(:));   

end%of function

function y = compute_Ax(x, p)
    x_f = psf2otf(x, p.img_size);
    y = otf2psf(p.m .* x_f, p.psf_size);
    y = y + p.lambda * x;
end
