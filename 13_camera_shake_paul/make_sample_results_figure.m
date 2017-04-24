imNums = [2 4];
kerNums = [7 4];
fontSize = 12;
figure(1); clf(1);

for ni = 1:length(imNums)
   imNum = imNums(ni);
   kerNum = kerNums(ni);
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

maxc = 0.1;
scl = @(x) max(0,x).^0.6;
ab = @(x) addborder(scl(x),2,maxc,'outer');
ps = {'SV',0.002,'SH',0.002,'ML',0.01,'MR',0.01,'MB',0.00,'MT',0.08};
 xLeft = 0; yTop = 0; xSize = 7; ySize = 2.8;
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

numRows = length(imNums);
numCols = 5;

subaxist(numRows,numCols,ni,1,ps);
   imagesc(x);
   if ni == 1, title('true','FontSize',fontSize); end
   colormap('gray'); axis image;
   set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
subaxist(numRows,numCols,ni,2,ps);
   imagesc(y);
   colormap('gray'); axis image;
   set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
   if ni == 1, title('blurred','FontSize',fontSize); end
   make_inset_axes;
   %addborder(img1, t, c, stroke)
   imagesc(ab(k)); axis image;
   set(gca,'XTick',[],'YTick',[],'CLim',[0 maxc]);
subaxist(numRows,numCols,ni,3,ps);
   imagesc(B_est_images{imNum,kerNum});
   colormap('gray'); axis image;
   set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
   if ni == 1, title('Babacan','FontSize',fontSize); end
   make_inset_axes;
   imagesc(ab(est_kernels{imNum,kerNum})); axis image;
   set(gca,'XTick',[],'YTick',[],'CLim',[0 maxc]);
subaxist(numRows,numCols,ni,4,ps);
   imagesc(fe.ex); colormap('gray'); axis image;
   set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
   if ni == 1, title('Levin','FontSize',fontSize); end
   make_inset_axes;
   imagesc(ab(flipud(fliplr(fe.k))));
   set(gca,'XTick',[],'YTick',[],'CLim',[0 maxc]); axis image;
subaxist(numRows,numCols,ni,5,ps);
   imagesc(ours.xEst); colormap('gray'); axis image;
   set(gca,'CLim',[0 1],'XTick',[],'YTick',[]);
   if ni == 1, title('Ours','FontSize',fontSize); end
   make_inset_axes;
   imagesc(ab(ours.kEst)); axis image;
   set(gca,'XTick',[],'YTick',[],'CLim',[0 maxc]);


end

%{
%ps = {'SV',0.01,'SH',0.01,'ML',0.01,'MR',0.01,'MB',0.05,'MT',0.05};%,...
   %'PL',0.02,'PR',0.02,'PB',0.02,'PT',0.02};
subaxist(numRows,numCols,1,6,ps);
%set(gca,'Position',[0.87 0.3 0.5 0.6]);
   err_bins = 10:10:150;
   nb = histc(ssdes_b(:),err_bins);
   ns = histc(errs(:),err_bins);
   nf = histc(feErrs(:),err_bins);
   hold on;
   plot(err_bins,100*cumsum(nb)/32,'b--','LineWidth',2);
   plot(err_bins,100*cumsum(nf)/32,'g-.','LineWidth',2);
   plot(err_bins,100*cumsum(ns)/32,'r-','LineWidth',2);

   xlabel('error (SSE)');
   ylabel('% below SSE');
   %plot(err_bins,[cumsum(nb) cumsum(nf) cumsum(ns) ]/32,'LineWidth',2);
   %legend({'Babacan','Levin','Ours'},'FontSize',fontSize,'Location','SouthOutside');
   set(gca,'YLim',100*[0 1.2],'YTick',100*(0.2:0.2:1));
subaxist(numRows,numCols,2,6,ps);
% timing graph: min, max, median
%}
print(1,'~/Desktop/2012_12_camL0/test.eps','-depsc','-r300');
