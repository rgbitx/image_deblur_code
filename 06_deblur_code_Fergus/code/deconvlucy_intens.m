function out = deconvlucy_intens(in,k,its)
  
  if (size(in,3)==1)
    error('Color images only');
  end
  
  % YIQ conversion
  q = rgb2ntsc(in);
  
  % deblur only Y channel
  q2 = deconvlucy(q(:,:,1),k,its);
  
  % make new YIQ image
  q3 = cat(3,q2,q(:,:,2),q(:,:,3));
  
  % back to RGB
  out = ntsc2rgb(q3);
  
  
