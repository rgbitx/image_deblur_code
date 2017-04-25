function out = deblurring_lr_whole_with_edge_without_lr(y, kernel, lr_wei, structure, thresh)
%----------------------------------------------------
% Deblurring with Low Rank method
% Date: April 14th, 2013
% Author: Jinshan Pan, jspan@mail.dlut.edu.cn
%----------------------------------------------------
%% for edges
dxf = [-1,1;0,0];
dyf = [-1,0;1,0];
% Sx = imfilter(structure,dxf, 'conv', 'circular');
% Sy = imfilter(structure,dyf, 'conv','circular');
Sx =  fft2(dxf, size(structure,1), size(structure,2));
Sx = ifft2(Sx.*fft2(structure));
Sy =  fft2(dyf, size(structure,1), size(structure,2));
Sy = ifft2(Sy.*fft2(structure));
S_mag = sqrt(Sx.^2+Sy.^2);
Sx = Sx.*heaviside_function(S_mag,thresh);
Sy = Sy.*heaviside_function(S_mag,thresh);
Sx(1:2,:) =0;Sx(end-1:end,:) = 0;Sx(:,1:2) =0;Sx(:,end-1:end) = 0;
Sy(1:2,:) =0;Sy(end-1:end,:) = 0;Sy(:,1:2) =0;Sy(:,end-1:end) = 0;
%% parameter setting
sigma = 25;
wei_grad = 5e-5;
%% End parameter setting
y = y-min(y(:));
y = y./max(y(:));
y_lr= y.*255;
%x_lr=Image_LASSC_Denoising(y_lr,y_lr,sigma);
x_lr = zeros(size(y));
%% Restoration
%  lr_wei = 0.01;
 kernel = flp(kernel);
 x_lr = double(x_lr)/255;
 out = deblurring_lr_with_edges(y, kernel, x_lr, lr_wei, Sx, Sy, wei_grad);
%  %%
%   for ii = 1:2
%      y_lr= out.*255;
%      x_lr=Image_LASSC_Denoising(y_lr,y_lr,sigma);
%      x_lr = double(x_lr)/255;
%      out = deblurring_lr_with_edges(y, kernel, x_lr, lr_wei, Sx, Sy, wei_grad);
%  end
%  




