function S = create_greenspan_settings(varargin)

% Author: Bryan Russell
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

% CREATE_GREENSPAN_SETTINGS - Creates a data structure containing the
% various parameter settings for the Greenspan nonlinear enhancement
% algorithm.  
%
% S = CREATE_GREENSPAN_SETTINGS creates the default settings and stores
% them in S.
%
% S = CREATE_GREENSPAN_SETTINGS('P1',SET1,'P2',SET2,...) sets the
% parameters P1 to be SET1, etc.  The settings that can be changed are as
% follows:
%
%   'lo_filt' - Lowpass filter, which can be any kernel or one of 'gauss'
%               or 'binomial5' (default 'binomial5')
%   'c' - Clipping constant (default 0.4).  You may also specify
%         'derived', which sets it to be the derived value, 0.45.  This
%         constant must be in [0,1].  This constant provides the
%         nonlinearity and augments the frequency content of the image.
%   's' - Scaling constant (default 5).  You may also specify 'derived'
%         which sets it to be the derived vaule, 3.  This constant
%         provides the sharp slope, thus reducing the blurring effect,
%         yet augmenting the ringing side-effect.
%   'bp' - Bandpass filtering flag, which can be either 0 or 1 (default
%          1). 
%   'factor' - Zoom by a factor of 2^factor (default is 1).
%
% See also GREENSPAN.
%
% <brussell@csail.mit.edu>
  
  S.lo_filt = binomialFilter(5);
  S.lo_filt = S.lo_filt*S.lo_filt';
  S.c = 0.4;
  S.s = 5;
  S.bp = 1;
  S.factor = 1;
  S = parse_args(S,varargin{:});
  
function S = parse_args(S,varargin)
  n = length(varargin);
  
  for k = 1:2:n-1
    param = varargin{k};
    setting = varargin{k+1};
    
    switch lower(param)
     case 'lo_filt'
      if ischar(setting)
	switch lower(setting)
	 case 'gauss'
	  S.lo_filt = fspecial('gauss',[3 3],1);
	 case 'binomial5'
	  S.lo_filt = binomialFilter(5);
	  S.lo_filt = S.lo_filt*S.lo_filt';
	 otherwise
	  error('Invalid setting of lo_filt');
	end
      else
	S.lo_filt = setting;
      end
     case 'c'
      if ischar(setting) & strcmp(lower(setting),'derived')
	S.c = 0.45;
      else
	S.c = setting;
      end
     case 's'
      if ischar(setting) & strcmp(lower(setting),'derived')
	S.s = 3;
      else
	S.s = setting;
      end
     case 'bp'
      S.bp = setting;
     case 'factor'
      S.factor = setting;
     otherwise
      error('Invalid parameter');
    end
  end
