
% opts.kernel_size = 35;
opts.kernel_size = 25;
% set kernel_est_win to be the window used for estimating the kEst - if
% this option is empty, whole image will be used
% opts.kernel_est_win = [428,576,535,638];
opts.kernel_est_win=[];
% set initial downsampling size for really large images
opts.prescale = 1;

% set this to 1 for no gamma correction - default 1.0
opts.gamma_correct = 1.0;

% kernel initialization at coarsest level
% 0 = uniform; 1 = vertical bar; 2 = horizontal bar; 3 = tiny 2-pixel
% wide kernel at coarsest level
opts.kernel_init = 3;
opts.medFilter = false;

% non-blind settings
opts.nb_lambda = 3000;
opts.nb_alpha = 1.0;
opts.use_ycbcr = 0;
fn = [];

% dataSource = 'levin_benchmark'; % 'field' or 'levin_benchmark' images
dataSource = 'field'
%%
switch dataSource
   case 'field'
      % parameter settings for real images:
      % real images conform less well to the uniform blur assumption,
      % but they are also larger and more sparse. we prefer a "quick
      % and dirty" approach with sparser initialization, less burn-in,
      % and a faster sparsity relaxation schedule.
%       fn = 'pietro.tif';
      fn = 'image3.png';
      %fn = 'mukta.jpg';
      %fn = 'fishes.jpg';
      %fn = 'lyndsey.tif';
      %  opts.medFilter = true; % remove lyndsey's salt-pepper noise
      %  opts.kernel_est_win = [135 55 1170 812]; % for lyndsey, to speed up computation

      opts.xk_iter = 30;     % x-k iterations in initial multiscale stages
      opts.final_scale_xk_iter = 60; % x-k iterations in final stage
      opts.burnIn = 10;      % x-k iterations in initial stage
      opts.stageInd = 5;     % iterations per increase in sparsity
      opts.pa_target = 1;    % pred/actual L1/L2 target
      opts.tFac0 = 0.10;     % sparsity factor for initial stage
      opts.tJumpFac = 1.10;  % jump in sparsity factor per stage
      opts.x_in_iter = 6;    % x inner iterations
      opts.k_in_iter = 6;    % k inner iterations
      opts.showFigs = true;  % show intermediate kernels, edge maps, residuals
      opts.decon_kernel = true; % use internal bregman deconv in ms_blind_deconv
      
      tic;
      [kEst,yorig,deblur] = ms_blind_deconv(fn, opts);
      etime = toc
      %%
      figure(26);
      ax(1)=subplot(121);
      imagesc(uint8(255*yorig));
      set(gca,'CLim',[0 255]); axis image;
      title('original');
      ax(2)=subplot(122);
      imagesc(uint8(255*deblur));
      set(gca,'CLim',[0 255]); axis image;
      title('deblurred');
      linkaxes(ax);
      figure(27);
      imagesc(kEst);
      
   case 'levin_benchmark'
      % conservative settings for levin set
      opts.xk_iter = 30;
      opts.final_scale_xk_iter = 180;
      opts.burnIn = 30;      % x-k iterations in initial stage
      opts.stageInd = 10;    % iterations per increase in sparsity
      opts.pa_target = 1;    % pred/actual L1/L2 target
      opts.tFac0 = 0.15;     % sparsity factor for initial stage
      opts.tJumpFac = 1.10;  % jump in sparsity factor per stage
      opts.x_in_iter = 1;    % x inner iterations
      opts.k_in_iter = 6;    % k inner iterations
      opts.showFigs = false;  % show intermediate kernels, edge maps, residuals
      opts.decon_kernel = false; % use internal bregman deconv in ms_blind_deconv
      
      imNum = 2; kerNum = 6;
      
      % load image and kernel
      imDir = './data/';
      imFile = ['im0' num2str(imNum) '_ker0' num2str(kerNum) '.mat'];
      imPath = [imDir imFile];
      load(imPath);
      kTrue = rot90(f,2);
      
      % show image
      figure(26);
      colormap('gray');
      subplot(131);
      imagesc(x); axis image;
      subplot(132);
      imagesc(kTrue); axis image;
      subplot(133);
      imagesc(y); axis image;
      opts.blur = y;
      
      % kernel estimation
      tic;
      kEst = ms_blind_deconv(fn, opts);
      etime = toc
      
      % deconvolve with estimated (kEst) and ground truth (kTrue)
      % using Levin's deconvolution
      xEst = deconvSps(y,rot90(kEst,2),0.0068,70);
      x0Est = deconvSps(y,rot90(kTrue,2),0.0068,70);
      
      % registration and error measurement:
      % find a translation that best matches up xEst and x,
      % get the sum of squared deviations between the matched up images
      err0 = comp_upto_shift(x0Est,x)
      errEst = comp_upto_shift(xEst,x)
      errRat = errEst/err0
      
      xMax = max(x(:));
      kMax = 0.25;
      figure(1000); colormap('gray');
      subplot(231);
      imagesc(x); set(gca,'CLim',[0 xMax]); axis image;
      title('true');
      subplot(232);
      imagesc(x0Est); set(gca,'CLim',[0 xMax]); axis image;
      title('decon true kEst');
      subplot(233);
      imagesc(xEst); set(gca,'CLim',[0 xMax]); axis image;
      title('decon est kEst');
      subplot(234);
      imagesc(y); set(gca,'CLim',[0 xMax]); axis image;
      title('blurred');
      p = 0.7;
      subplot(235);
      imagesc(kTrue.^p); set(gca,'CLim',[0 kMax]); axis image;
      subplot(236);
      imagesc(kEst.^p); set(gca,'CLim',[0 kMax]); axis image;
end
