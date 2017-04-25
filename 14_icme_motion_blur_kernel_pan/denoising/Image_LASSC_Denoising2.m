%----------------------------------------------------
% Blur: 9x9 Uniform kernel; AGWN std Varance = 2^0.5
% Data: May 20th, 2010
% Author: Weisheng Dong, wsdong@mail.xidian.edu.cn
%----------------------------------------------------
function  [xx,r]=Image_LASSC_Denoising2(y, x,  nSig, win, iter)
pre           =   'LASSC_';
par.nSig      =   nSig;
par.win=win;
par.K=iter;
if nSig<=15
    %par.win       =   6;
    par.nblk      =   40;
    par.c1        =   2.8*sqrt(2);   % 1.7
    par.lamada    =   0.63;
    par.w         =   0.23;
    %par.K         =   10;
elseif nSig <= 30
    %par.win       =   7;
    par.nblk      =   60;
    par.c1        =   2.9*sqrt(2);   % 2.6
    par.lamada    =   0.65;
    par.w         =   0.23;
    %par.K         =   11;
elseif nSig<=50
    %par.win       =   8;
    par.nblk      =   75;
    par.c1        =   3.0*sqrt(2);   % 2.6
    par.lamada    =   0.67;
    par.w         =   0.23;    
    %par.K         =   12;
else
    %par.win       =   9;
    par.nblk      =   90;
    par.c1        =   3.1*sqrt(2);   % 1.6
    par.lamada    =   0.64;
    par.w         =   0.23;    
    %par.K         =   14;    
end
par.step      =   min(4, par.win-1);

    
    par.I        =   [x fliplr(x);flipud(x) flipud(fliplr(x))];
    par.nim      =   [y fliplr(y);flipud(y) flipud(fliplr(y))];
    
    [xx_buf PSNR SSIM]   =   LASSC_Denoising( par );
    xx=(xx_buf(1:end/2,1:end/2)+fliplr(xx_buf(1:end/2,end/2+1:end))...
        +flipud(xx_buf(end/2+1:end,1:end/2))+fliplr(flipud(xx_buf(end/2+1:end,end/2+1:end))))/4;
    
   %disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n', im_dir(i).name, PSNR, SSIM) );
    

