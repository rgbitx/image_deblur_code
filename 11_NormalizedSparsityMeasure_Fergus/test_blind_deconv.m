clear
close all
opts.kernel_size = 25;

% set kernel_est_win to be the window used for estimating the kernel - if
% this option is empty, whole image will be used
opts.kernel_est_win = []; 

% set initial downsampling size for really large images
opts.prescale = 1;

% This is the weight on the likelihood term - it should be decreased for
% noisier images; decreasing it usually makes the kernel "fatter";
% increasing makes the kernel "thinner".
opts.min_lambda = 100;

% Kernel regularization weight
opts.k_reg_wt = 0.0;

% set this to 1 for no gamma correction - default 1.0
opts.gamma_correct = 1.0;

% threshold on fine scale kernel elements
opts.k_thresh = 0.0;

% kernel initialiazation at coarsest level
% 0 = uniform; 1 = vertical bar; 2 = horizontal bar; 3 = tiny 2-pixel
% wide kernel at coarsest level
opts.kernel_init = 3; 

% delta step size for ISTA updates; increasing this delta size is not a
% good idea since it may cause divergence. On the other hand decreasing
% it too much will make convergence much slower. 
opts.delta = 0.001;

% inner/outer iterations for x estimation
opts.x_in_iter = 2; 
opts.x_out_iter = 2;

% maximum number of x/k alternations per level; this is a trade-off
% between performance and quality.
opts.xk_iter = 21;

% non-blind settings
opts.nb_lambda = 3000;
opts.nb_alpha = 1.0;
opts.use_ycbcr = 1;

% Note that the min_lambda parameter should be varied for different images to
% give better results.
fn = 'lyndsey.tif'; opts.kernel_est_win = [335 275 1170 712]; opts.min_lambda = 100;
% fn = 'pietro.tif'; opts.kernel_est_win = [141 123 601 828]; opts.min_lambda = 100;
% fn = 'mukta.jpg'; opts.kernel_size = 27; opts.use_ycbcr = 0; opts.min_lambda = 200;

% From Cho/Lee et. al. SIGGRAPH Asia 2009
% fn = 'image5.png'; opts.kernel_size = 35; opts.min_lambda = 200;

[blur, deblur, kernel, opts] = ms_blind_deconv(fn, opts);
figure
imshow(blur)
figure
imshow(deblur)




