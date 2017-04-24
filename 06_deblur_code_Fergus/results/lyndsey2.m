
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User specified parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User selected region (xmin xmax ymin ymax)
AXIS = [121   315   681   874];
% max size of blur kernel
BLUR_KERNEL_SIZE = 25;
% inital value of kernel
%FIRST_INIT_MODE_BLUR = 'hbar';
FIRST_INIT_MODE_BLUR = 'vbar';
%FIRST_INIT_MODE_BLUR = 'delta';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image to deblur
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename of image to deblur
obs_im = imread('../images/lyndsey2.jpg');  

% downsample image before we start to keep image size managable
PRESCALE = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Richardson - Lucy related parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of lucy richardson iterations to use
% default = 10 but more for long thin blurs...
LUCY_ITS = 10;

% 0=use full res kernel, 1=use 1 scale coarser etc.
SCALE_OFFSET = 0; 

% threshold on kernel intensities (helps to remove noise from kernel)
% is % of max value in kernel
KERNEL_THRESHOLD = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other semi-important parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 = use rgb2gray; 1 = just use red channel, 2 = green; 3 = blue
SPECIFIC_COLOR_CHANNEL = 0; 

% type of prior to use
%PRIOR_TYPE = 'whiteboard';
PRIOR_TYPE = 'street';

% intensity value above which counts as saturation
SATURATION_THRESHOLD = 250;
% camera type 
CAMERA_TYPE = 'unknown';

% gamma correction (set to 1 for RAW images)
GAMMA_CORRECTION = 2.2;

% use automatic patch selector
AUTOMATIC_PATCH = 0;
% run synthetic experiment or not
SYNTHETIC = 0;
% apply blur kernel priors or not
BLUR_LOCK = 1;
% re-center kernel after every scale
CENTER_BLUR = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters - dont alter unless you really know what you are doing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RESIZE_STEP = sqrt(2); %% if ~=2 must use matlab..
NUM_SCALES = ceil(-log(3/BLUR_KERNEL_SIZE)/log(sqrt(2)))+1;

if ~AUTOMATIC_PATCH
  PATCH_SIZE = [AXIS(2)-AXIS(1) AXIS(4)-AXIS(3)];
  PATCH_LOCATION = [AXIS(1) AXIS(3)];
else
  PATCH_SIZE = [256 256];
end

BLUR_MASK = 0;
BLUR_MASK_VARIANCES = [50 0.02 5 0.2];

GRADIENT_MODE = 'haar';
%GRADIENT_MODE = 'steer';

%RESIZE_MODE = 'matlab_nearest';
RESIZE_MODE = 'matlab_bilinear';
%RESIZE_MODE = 'matlab_bicubic';
%RESIZE_MODE = 'binom5';

EXTRA_BLUR = 0;
EXTRA_BLUR_SIZE = 5;
EXTRA_BLUR_SIGMA = 0.5;

BLUR_MASK = 0;
BLUR_MASK_VARIANCES = [100 0.01];

AUTOMATIC_PATCH_CENTER_WEIGHT = 5; 
SATURATION_MASK = 1;
FFT_MODE = 1;

if (GAMMA_CORRECTION==1)
  INTENSITY_SCALING = 1/256;
else
  INTENSITY_SCALING = 1; 
end

%UPSAMPLE_MODE = 'bill_filter';
%UPSAMPLE_MODE = 'greenspan';
UPSAMPLE_MODE = 'matlab_bilinear';

RESCALE_THEN_GRAD = 0; %% order to do rescaling/gradient operations (only
                       %for non-eero functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall system parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEBUG = 0;
SYNTHETIC = 0;
COLOR = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization parameters to get x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIRST_INIT_MODE_BLUR = 'lucy';
%FIRST_INIT_MODE_BLUR = 'reg';
%FIRST_INIT_MODE_BLUR = 'variational';
%FIRST_INIT_MODE_BLUR = 'direct';
%FIRST_INIT_MODE_BLUR = 'true';
%FIRST_INIT_MODE_BLUR = 'random';
%FIRST_INIT_MODE_BLUR = 'updown';
%FIRST_INIT_MODE_BLUR = 'vbar';
%FIRST_INIT_MODE_BLUR = 'delta';

%FIRST_INIT_MODE_IMAGE = 'lucy';
%FIRST_INIT_MODE_IMAGE = 'reg';
FIRST_INIT_MODE_IMAGE = 'variational';
%FIRST_INIT_MODE_IMAGE = 'direct';
%FIRST_INIT_MODE_IMAGE = 'true';
%FIRST_INIT_MODE_IMAGE = 'random';
%FIRST_INIT_MODE_IMAGE = 'slight_blur_obs';
%FIRST_INIT_MODE_IMAGE = 'updown';
%FIRST_INIT_MODE_IMAGE = 'nearest';
%FIRST_INIT_MODE_IMAGE = 'greenspan';

%INIT_MODE_BLUR = 'lucy';
%INIT_MODE_BLUR = 'reg';
%INIT_MODE_BLUR = 'variational';
INIT_MODE_BLUR = 'direct';
%INIT_MODE_BLUR = 'true';
%INIT_MODE_BLUR = 'updown';

%INIT_MODE_IMAGE = 'lucy';
%INIT_MODE_IMAGE = 'reg';
%INIT_MODE_IMAGE = 'variational';
INIT_MODE_IMAGE = 'direct';
%INIT_MODE_IMAGE = 'true';
%INIT_MODE_IMAGE = 'updown';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLUR_PRIOR = 1; % Exponential
IMAGE_PRIOR = 0; % Gaussian

if BLUR_LOCK
  BLUR_COMPONENTS = 4;
else
  BLUR_COMPONENTS = 3;
end

IMAGE_COMPONENTS = 4;
BLUR_UPDATE_FREQ = 1;
INIT_PRESCISION = 1e4;
NOISE_INIT = 1;
CONVERGENCE = 5e-4;
MAX_ITERATIONS = 50000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post iteration operations on blur kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%POST_PROCESS_BLUR = 'noise_removal';
POST_PROCESS_BLUR = 'none';
POST_PROCESS_THRESHOLD = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstsruct image mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMAGE_RECONSTRUCTION = 'lucy';
%IMAGE_RECONSTRUCTION = 'reg';
%IMAGE_RECONSTRUCTION = 'variational';


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior file name root
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CAMERA_TYPE,'unknown')
  if strcmp(PRIOR_TYPE,'street')
    prior_name = ['../priors/linear_street_',num2str(IMAGE_COMPONENTS)];
  elseif strcmp(PRIOR_TYPE,'whiteboard')
    prior_name = ['../priors/linear_whiteboard_',num2str(IMAGE_COMPONENTS)];
  else
    error('Unknown prior type');
  end
else
  prior_name = ['../priors/',CAMERA_TYPE,'_street_',num2str(IMAGE_COMPONENTS)];
end
load(prior_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel for synthetic mode or blur guess for real mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SYNTHETIC
  error('Please provide your own synthetic kernel - not included in release version of code');
  tmp = load('blur_kernel17.mat');
  blur_kernel = tmp.blur_kernel; % may need to change variable name
  blur_kernel = blur_kernel/sum(blur_kernel(:));
else
  blur_kernel = delta_kernel(BLUR_KERNEL_SIZE);
end

% end of options section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save all to current filename
fname = mfilename;
save(fname);

