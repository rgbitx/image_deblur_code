kesize = 35;
fontSize = 14;

dataDir = './data/';

% babacan results -- all in one file
load([dataDir 'babacan_results.mat']);

% shearer and levin results -- in 32 files
for imNum = 1:4
   for kerNum = 1:8
      imFile = ['im0' num2str(imNum) '_ker0' num2str(kerNum) '.mat'];
      imPath = [dataDir imFile];
      load(imPath,'f','x','y');
      k = fliplr(flipud(f));
      % shearer results
      deconFile = ['isep_' imFile];
      deconPath = [dataDir deconFile];
      ours = load(deconPath,'xEst','errRat','kEst','err0','etime');
      k_size = size(k);
      kCent = floor(size(ours.kEst)/2)+1;
      [X,Y] = meshgrid(1:size(ours.kEst,1));
      kCent(1)=round(sum(vec(ours.kEst.*Y))); 
      kCent(2)=round(sum(vec(ours.kEst.*X))); 
      kRad = floor(k_size(1)/2);
      ours.kEst = ours.kEst(...
         kCent(1)-kRad:kCent(1)+kRad,...
         kCent(2)-kRad:kCent(2)+kRad);
         
      % levin results
      fileEnd = ['_im' num2str(imNum) ...
         '_ker' num2str(kerNum) '.mat'];
      feFileStem = 'diagfe_filt_sps';
      feFile = [feFileStem fileEnd];      
      fePath = [dataDir feFile];
      fe = load(fePath,'ex','k','ssde');

      ps = {'SV',0.03,'SH',0.03,'ML',0.05,'MR',0.05,'MB',0.05,'MT',0.05};
       xLeft = 0; yTop = 0; xSize = 10; ySize = 2;
      %figure(figNum); close(figNum); figure(figNum);
      set(gcf,'PaperUnits','inches');
      set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
      ppos = get(gcf,'PaperPosition');

      su = get(gcf,'Units');
      pu = get(gcf,'PaperUnits'); 
      set(gcf,'Units',pu);
      spos = get(gcf,'Position');
      set(gcf,'Position',[spos(1) spos(2) [ppos(3) ppos(4)]*1.22523994]);%*1.22523994)
      set(gcf,'Units',su) 

      numRows = 1;
      ab = @(x) addborder(x,3,max(x(:)),'outer');
      figure(1);
      subaxist(numRows,5,1,1);
      imagesc(x); 
         colormap('gray'); axis image; 
         set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
      subaxist(numRows,5,1,2);
         imagesc(y); 
         colormap('gray'); axis image;
         set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
         make_inset_axes;
         %addborder(img1, t, c, stroke)
         imagesc(ab(k));
         set(gca,'XTick',[],'YTick',[]);
      subaxist(numRows,5,1,3);
         imagesc(ours.xEst); colormap('gray'); axis image;
         set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
         title('Ours');
         make_inset_axes;
         imagesc(ab(ours.kEst));
         set(gca,'XTick',[],'YTick',[]);
      subaxist(numRows,5,1,4);
         imagesc(B_est_images{imNum,kerNum}); 
         colormap('gray'); axis image;
         set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
         title('Babacan');
         make_inset_axes;
         imagesc(ab(est_kernels{imNum,kerNum}));
         set(gca,'XTick',[],'YTick',[]);
      subaxist(numRows,5,1,5);
         imagesc(fe.ex); colormap('gray'); axis image;
         set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
         title('Levin');
         make_inset_axes;
         imagesc(ab(flipud(fliplr(fe.k))));
         set(gca,'XTick',[],'YTick',[]);
      pause;
      
   end
end
