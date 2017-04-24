function [mx_new,me_new]=move_level(mx,me,K,L,M,N,mode,resize_step,center)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology


  %%% Get number of image channels
  mxc = size(mx,3);
  
  if (nargin==8)
    center = 0; %% no centering of kernel by default
  end

  if (center)

    fprintf('Centering kernel\n');
    me = me / sum(me(:));
    %% get centre of mass
    mu_y = sum([1:size(me,1)] .* sum(me,2)');
    mu_x = sum([1:size(me,2)] .* sum(me,1));    
  
    %% get mean offset
    offset_x = round( floor(size(me,2)/2)+1 - mu_x );
    offset_y = round( floor(size(me,1)/2)+1 - mu_y );
    
    %% make kernel to do translation
    shift_kernel = zeros(abs(offset_y*2)+1,abs(offset_x*2)+1);
    shift_kernel(abs(offset_y)+1+offset_y,abs(offset_x)+1+offset_x) = 1;
    
    %% shift both image and blur kernel
    me = conv2(me,shift_kernel,'same');
    
    for c=1:mxc
      mx(:,:,c) = conv2(mx(:,:,c),flipud(fliplr(shift_kernel)),'same');
    end
    
  end
  
if any(findstr(mode,'matlab')) 
  
  if strcmp(mode,'matlab_nearest')
    mx_new = imresize(mx,[M N],'nearest');
    me_new = imresize(me,[K L],'nearest');    
  elseif strcmp(mode,'matlab_bilinear')
    mx_new = imresize(mx,[M N],'bilinear');
    me_new = imresize(me,[K L],'bilinear');   
  elseif strcmp(mode,'matlab_bicubic')
    mx_new = imresize(mx,[M N],'bicubic');
    me_new = imresize(me,[K L],'bicubic');   
  else
    mx_new = imresize(mx,[M N],'bilinear');
    me_new = imresize(me,[K L],'bilinear');
  end
elseif strcmp(mode,'greenspan')
    
    s = create_greenspan_settings;
    
    [mx_new,lo] = greenspan(mx,s);
    
    if (resize_step~=2)
      mx_new = imresize(mx_new,[M N],'bilinear');
    end

     me_new = imresize(me,[K L],'bilinear');

elseif strcmp(mode,'bill_filter')
  
  % upsample as per usual
  mx_new = imresize(mx,[M N],'bilinear');
  me_new = imresize(me,[K L],'bilinear');
 
    
  %% make bill's filter
  bill_coeffs = [-0.0625 -0.25 1.625 -0.25 -0.0625];
  bill_filter = bill_coeffs' * bill_coeffs;
  %% convolve with image
  for c=1:mxc
    mx_new(:,:,c) = conv2(mx_new(:,:,c),bill_filter,'same');
  end
  
else
  
  %%%% NOT SURE IF THIS WORKS FOR COLOR
  
  if (mxc>1) %% temporary
    error foo
  end
  
  mx_new = upBlur(mx,1,mode);
  me_new = upBlur(me,1,mode);
    
  %%% check for size
  me_new = me_new(1:K,1:L);
  mx_new = mx_new(1:M,1:N);
  
end  




%%% ensure blur kernel is normalized
mfactor = sum(me_new(:));
me_new = me_new / mfactor;

%%% now 
%mx_new = mx_new * mfactor;
