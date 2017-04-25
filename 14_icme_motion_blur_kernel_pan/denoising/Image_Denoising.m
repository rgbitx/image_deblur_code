%----------------------------------------------------
% Blur: 9x9 Uniform kernel; AGWN std Varance = 2^0.5
% Data: May 20th, 2010
% Author: Weisheng Dong, wsdong@mail.xidian.edu.cn
%----------------------------------------------------
function  Image_Denoising( nSig, Out_dir, In_dir )
pre           =   'LASSC_';
par.nSig      =   nSig;
if nSig<=15
    par.win       =   6;
    par.nblk      =   40;
    par.c1        =   2.8*sqrt(2);   % 1.7
    par.lamada    =   0.63;
    par.w         =   0.23;
    par.K         =   10;
elseif nSig <= 30
    par.win       =   7;
    par.nblk      =   60;
    par.c1        =   2.9*sqrt(2);   % 2.6
    par.lamada    =   0.65;
    par.w         =   0.23;
    par.K         =   11;
elseif nSig<=50
    par.win       =   8;
    par.nblk      =   75;
    par.c1        =   3.0*sqrt(2);   % 2.6
    par.lamada    =   0.67;
    par.w         =   0.23;    
    par.K         =   12;
else
    par.win       =   9;
    par.nblk      =   90;
    par.c1        =   3.1*sqrt(2);   % 1.6
    par.lamada    =   0.64;
    par.w         =   0.23;    
    par.K         =   14;    
end
par.step      =   min(6, par.win-1);
fpath         =   fullfile(In_dir, '*.tif');
im_dir        =   dir(fpath);
im_num        =   length(im_dir);
cnt           =   0;
sum_psnr      =   0;
sum_ssim      =   0;
time0         =   clock;
fn_txt        =   strcat( pre, 'PSNR_SSIM.txt' ); 
fd_txt        =   fopen( fullfile(Out_dir, fn_txt), 'wt');

for i = 1:im_num 
    
    par.I        =   double( imread(fullfile(In_dir, im_dir(i).name)) );
    par.nim      =   par.I + nSig*Gen_noise(In_dir, im_dir, i);
    
    [im PSNR SSIM]   =   LASSC_Denoising( par );
    sum_psnr    =  sum_psnr + PSNR;
    sum_ssim    =  sum_ssim + SSIM;    
    
    fname            =   strcat(pre, im_dir(i).name);
    imwrite(im./255, fullfile(Out_dir, fname));
    disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n', im_dir(i).name, PSNR, SSIM) );
    fprintf(fd_txt, '%s :  PSNR = %2.2f  SSIM = %2.4f\n', im_dir(i).name, PSNR, SSIM);
    cnt   =  cnt + 1;
end
fprintf(fd_txt, '\n\nAverage :  PSNR = %2.2f  SSIM = %2.4f\n', sum_psnr/cnt, sum_ssim/cnt);
fclose(fd_txt);
disp(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));
return;


function nim  =  Gen_noise( In_dir, im_dir, i )
randn('seed',0);
for ii=1:i
    im        =   imread(fullfile(In_dir, im_dir(ii).name));
    nim       =   randn(size(im));
end
return;

