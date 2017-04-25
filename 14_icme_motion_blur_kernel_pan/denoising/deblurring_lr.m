function out = deblurring_lr(y, kernel, x_lr, lr_wei)
%----------------------------------------------------
% Deblurring with Low Rank method
% Date: April 14th, 2013
% Author: Jinshan Pan, jspan@mail.dlut.edu.cn
%----------------------------------------------------
if(max(x_lr(:))>128)
    x_lr = double(x_lr)/255;
end
%% Boundary processing
H = size(y,1);    W = size(y,2);
y_pad = wrap_boundary_liu(y, opt_fft_size([H W]+size(kernel)-1));
tmp = zeros(size(y_pad));
tmp(1:H,1:W) = x_lr;
x_lr = tmp;
clear tmp;
 [N,M,D] = size(x_lr);
 sizeI2D = [N,M];
 KER = psf2otf(kernel,sizeI2D);
 Den_KER = abs(KER).^2;
 Normin_lr = fft2(x_lr);
 Normin = conj(KER).*fft2(y_pad);
 Denormin = Den_KER + lr_wei;
 FI = (Normin+lr_wei*Normin_lr)./Denormin;
 out = real(ifft2(FI));
 out = out(1:H,1:W,:);
 %%
%  out(out<0) = 0;
%  out(out>1) = 1;
%   out = out-min(out(:));
% out = out/max(out(:));
 %figure; imshow(I)


