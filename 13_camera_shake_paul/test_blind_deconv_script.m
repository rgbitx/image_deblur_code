
opts.kernel_size = 35;

% set kernel_est_win to be the window used for estimating the kernel - if
% this option is empty, whole image will be used
opts.kernel_est_win = []; 

% set initial downsampling size for really large images
opts.prescale = 1;

% set this to 1 for no gamma correction - default 1.0
opts.gamma_correct = 1.0;

% kernel initialiazation at coarsest level
% 0 = uniform; 1 = vertical bar; 2 = horizontal bar; 3 = tiny 2-pixel
% wide kernel at coarsest level
opts.kernel_init = 3; 

% maximum number of x/k alternations per level; this is a trade-off
% between performance and quality.
opts.xk_iter = 30;
opts.final_scale_xk_iter = 180;

% settings
opts.burnIn = 30;      % x-k iterations in initial stage
opts.stageInd = 10;    % iterations per increase in sparsity
opts.pa_target = 1;    % pred/actual L1/L2 target
opts.tFac0 = 0.15;     % sparsity factor for initial stage
opts.tJumpFac = 1.10;  % jump in sparsity factor per stage
opts.x_in_iter = 1;    % x inner iterations
opts.k_in_iter = 6;    % k inner iterations
opts.showFigs = false; % show intermediate kernels, edge maps, residuals

% use internal decon routine?
opts.decon_kernel = false;

% load data
imPath = [imDir imFile];
load(imPath);
k = rot90(f,2);
opts.blur = y;

%% kernel estimation

tic
kEst = ms_blind_deconv([], opts);
etime = toc;

%% compare to non-blind decon with ground truth PSF using Levin's decon alg

xEst = deconvSps(y,rot90(kEst,2),0.0068,70);
x0Est = deconvSps(y,rot90(k,2),0.0068,70);

err0 = comp_upto_shift(x0Est,x)
errEst = comp_upto_shift(xEst,x)

errRat = errEst/err0

deconFile = ['isep_' imFile];
deconPath = [deconDir deconFile];
save(deconPath,'err0','errEst','errRat','xEst','kEst','etime');
