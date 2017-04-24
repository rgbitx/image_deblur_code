function ExportFig(Fig_Num,fname,mode,Res)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper  
% Copyright 2006, Massachusetts Institute of Technology 

if (nargin==3)
   Res=150;
end

if (nargin==2)
   mode='eps';
   Res=150;
end


if strcmp(mode,'eps')
   
   ext_ind = findstr(fname,'.eps');
   
   if (~isempty(ext_ind))
      fname = fname(1:ext_ind-1);
   end
   cmd=['print(',num2str(Fig_Num),',''-depsc2'',','''-r',num2str(Res),''',''',fname,'.eps'')'];
   eval(cmd);

elseif (strcmp(mode,'jpg') | strcmp(mode,'jpeg'))

   ext_ind = findstr(fname,'.jp');
   
   if (~isempty(ext_ind))
      fname = fname(1:ext_ind-1);
   end
   cmd=['print(',num2str(Fig_Num),',''-djpeg'',''',fname,'.jpg'')'];
   
   %cmd=['print(',num2str(Fig_Num),',''-djpeg'',','''-r',num2str(Res),''',''',fname,'.jpg'')'];
   t=evalc(cmd);

elseif (strcmp(mode,'png'))

   ext_ind = findstr(fname,'.png');
   
   if (~isempty(ext_ind))
      fname = fname(1:ext_ind-1);
   end
   
   cmd=['print(',num2str(Fig_Num),',''-dpng'',','''-r',num2str(Res),''',''',fname,'.png'')'];
   evalc(cmd);

elseif (strcmp(mode,'tif') | strcmp(mode,'tiff'))

   ext_ind = findstr(fname,'.ti');
   
   if (~isempty(ext_ind))
      fname = fname(1:ext_ind-1);
   end
   
   cmd=['print(',num2str(Fig_Num),',''-dtiff'',','''-r',num2str(Res),''',''',fname,'.tif'')'];
   evalc(cmd);

else
   error('Unsupported export format');
end
