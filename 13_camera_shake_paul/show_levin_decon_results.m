kesize = 35;
fontSize = 14;

dataDir = './data/';

% babacan results -- all in one file
b=load([dataDir 'babacan_results.mat']);

% shearer and levin results -- in 32 files
for imNum = 1:4
   for kerNum = 1:8
      imFile = ['im0' num2str(imNum) '_ker0' num2str(kerNum) '.mat'];
      imPath = [dataDir imFile];
      load(imPath,'f');

      % shearer results
      deconFile = ['isep_' imFile];
      deconPath = [dataDir deconFile];
      load(deconPath,'errRat','kEst','err0','etime');
      
      % levin results
      fileEnd = ['_im' num2str(imNum) ...
         '_ker' num2str(kerNum) '.mat'];
      feFileStem = 'diagfe_filt_sps';
      feFile = [feFileStem fileEnd];      
      fePath = [dataDir feFile];
      
      errs(imNum,kerNum) = errRat*err0;
      load(fePath,'ssde');
      feErr = ssde;
      feErrs(imNum,kerNum) = feErr;
      
      errRats(imNum,kerNum) = errRat;
      feErrRats(imNum,kerNum) = feErr/err0;

      etimes(imNum,kerNum) = etime;
      
      k = rot90(f,2);
      disp(['imNum = ' num2str(imNum) '   kerNum = ' num2str(kerNum)' ...
         '   errRat = ' num2str(errRat)]);
      imInd1 = 1:kesize;
      imInds = 1+(imNum-1+1)*kesize : (imNum+1)*kesize;
      kerInds = 1+(kerNum-1)*kesize : kerNum*kesize;
      keBig(imInd1,kerInds) = rot90(centerpad(k,kesize*ones(1,2)),2);
      keBig(imInds,kerInds) = rot90(kEst,2);

      
   end
end
%%
p = 0.5;
figure(98);
imagesc(keBig.^p); colormap('gray'); axis image;
title('kernels from Shearer method');

figure(99); clf(99);
 xLeft = 0; yTop = 0; xSize = 3+3/8; ySize = 2.5;
%figure(figNum); close(figNum); figure(figNum);

ps = ...
    {'SV',0.06,'SH',0.05,...
     'ML',0.10,'MR',0.05,'MB',0.10,'MT',0.02,...
     'PL',0.01,'PR',0.01,'PT',0.01,'PB',0.01};
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
ppos = get(gcf,'PaperPosition');

su = get(gcf,'Units');
pu = get(gcf,'PaperUnits'); 
set(gcf,'Units',pu);
spos = get(gcf,'Position');
set(gcf,'Position',[spos(1) spos(2) ppos(3) ppos(4)]*1.22523994)
set(gcf,'Units',su) 


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
legend({'Babacan','Levin','Ours'},'FontSize',fontSize,'Location','SouthEast');
set(gca,'YLim',100*[0 1.2],'YTick',100*(0.2:0.2:1));
