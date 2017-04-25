clc;
clear;
close all;
addpath(genpath('image'));
addpath(genpath('Levin_deblurring_code'));
addpath(genpath('denoising'));
addpath(genpath('filter'));
addpath(genpath('TV_deconvolution_model'));
 %parameter setting
 lambda_kernel = 0.01;% 0.01
 lambda_smooth = 0.005;%
 lambda_texture = 1;
window_size = 3;


 %lambda_kernel_smooth = 0.001;% 0.001 to 0.1
%  prompt = {'Enter the smoothness weight of kernel:'};
%  dlg_title = 'Weight of kernel parameter';
%  num_lines = 1;
%  def = {'1e-3'};
%  lambda_kernel_smooth = inputdlg(prompt,dlg_title,num_lines,def);
%  lambda_kernel_smooth = str2double(lambda_kernel_smooth{1});
 lambda_kernel_smooth = 1e-4;
 
% prompt = {'Display every iterations:(Y/N):'};
% dlg_title = 'Input the size of kernel';
% num_lines = 1;
% def = {'Y'};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% if answer{1}=='Y'
%     display = 1;
% else
%     display = 0;
% end
display = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rmse_c = zeros(4,8);
rmse = zeros(4,8);
ssde_no_lr = zeros(4,8);
dirname = 'resdir_no_low_rank';
for i=1:4
    for j=[5,3,2,1,6,7,8,4]
        [i,j]
        i = 2
        j = 1
        eval(sprintf('load test_data/im%02d_ker%02d',i,j))
        mat_outname=sprintf('%s/im%d_ker%d.mat',dirname,i,j);
        d=dir(mat_outname);
        B_patch = y;
        ori_B =  y;
        [kernel_sizeh, kernel_sizew]= size(f);
        [our_deblur, deblur_levin, ssde, kernel] = deconv_main_multiple(B_patch, ori_B, lambda_kernel, lambda_smooth,...
            lambda_texture, window_size, kernel_sizeh, kernel_sizew, lambda_kernel_smooth, display,x);
        %%%
        kernel_c  = Cho_correct(kernel);
        kernel_c(kernel_c<0) = 0;
        kernel_c = kernel_c/sum(kernel_c(:));
        rmse_c(i,j) = sum((kernel_c(:) - f(:)).^2);
        rmse(i,j) = sum((kernel(:) - f(:)).^2);
        eval(sprintf('save  %s kernel deblur_levin our_deblur ssde' ,mat_outname));
        kk = kernel;
        kk = kk-min(kk(:));
        kk = kk/max(kk(:));
        sss = sprintf('im%d_ker%d_kernel.png',i,j);
        ssde_no_lr(i,j) = ssde;
        imwrite(kk,['levin_kernel\' sss]);
    end
end
save('resdir_no_low_rank\ssde_no_lr', 'ssde_no_lr');
    
    
    
    
    
    
    
    