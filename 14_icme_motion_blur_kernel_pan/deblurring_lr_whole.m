function out = deblurring_lr_whole(y, kernel, lr_wei)
%----------------------------------------------------
% Deblurring with Low Rank method
% Date: April 14th, 2013
% Author: Jinshan Pan, jspan@mail.dlut.edu.cn
%----------------------------------------------------
sigma = 25;
y = y-min(y(:));
y = y./max(y(:));
y_lr= y.*255;
x_lr=Image_LASSC_Denoising(y_lr,y_lr,sigma);
%% Restoration
%  lr_wei = 0.01;
 kernel = flp(kernel);
 x_lr = double(x_lr)/255;
 out = deblurring_lr(y, kernel, x_lr, lr_wei);
 %%
  for ii = 1:2
     y_lr= out.*255;
     x_lr=Image_LASSC_Denoising(y_lr,y_lr,sigma);
     x_lr = double(x_lr)/255;
     out = deblurring_lr(y, kernel, x_lr, lr_wei);
 end
 




