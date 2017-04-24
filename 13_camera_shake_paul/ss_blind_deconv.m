function [x, k] = ss_blind_deconv(y, x, k, opts)

   xk_iter     = opts.xk_iter;
   stageInd    = opts.stageInd;
   burnIn      = opts.burnIn;
   pa_target   = opts.pa_target;
   tFac0       = opts.tFac0;
   tJumpFac    = opts.tJumpFac;
   x_in_iter   = opts.x_in_iter;
   k_in_iter   = opts.k_in_iter;
   showFigs    = opts.showFigs;

   xk_iter = xk_iter+1;
   
   kLog = zeros(size(k,1),size(k,2),xk_iter);
   kLog(:,:,1) = k;
   
   timePerKernDisp = 5;
   timePerDataDisp = 5;
   t0k = tic;
   t0d = tic;
   
   L1g = @(x) sum(vec(normd([vec(x{1}), vec(x{2})],2,2)));
   L2g = @(x) norm([vec(x{1}); vec(x{2})],2);
   L12g = @(x) L1g(x)/L2g(x);
   
   khs = floor(size(k, 1) / 2);
   
   % arrays to hold costs
   lcost = zeros(xk_iter+1,1);
   pcost = zeros(xk_iter+1,1);
      
   ypcost = L12g(y)^2;
   
   sizes_x1 = [size(x{1}); size(x{2})];
   numels_x1 = [numel(x{1}); numel(x{2})];
   gc2vec2 = @(x) [x{1}(:),x{2}(:)];
   vec22gc = @(x) {reshape(x(:,1),sizes_x1(1,:)),...
      reshape(x(:,2),sizes_x1(2,:))};
   
   tFac = tFac0;
   t = min(numel(y{1}),max(4,tFac*ypcost));
   proj = @(x) vec22gc(vec_prox(@top_k_hard,gc2vec2(x),round(t)));
   conFun = @(x) sum(L2g(x) > 0);
   
   x = proj(x);
   res{1} = conv2(x{1}, k, 'same') - y{1};
   res{2} = conv2(x{2}, k, 'same') - y{2};
   
   iter = 1;
   
   lcost(iter) = (1/2)*(norm(res{1}(:))^2 + norm(res{2}(:))^2);
   pcost(iter) = conFun(x);
      
   fprintf('%5s %13s %13s %13s %13s\n', ...
      'iter','fit','nnz','x/y L1/L2','kx/y L1/L2');
   
   ti = 2;
   for iter = 2:xk_iter
      
      %%% x update: 
      %%% projected gradient descent with an optimal
      %%% step size and variable sparsity.      
      
      max_ls_iter = 4;
      x_in_iter = 1;
      end_x_iter = false;
      optTol = 0;
      for inn_iter=1 : x_in_iter
         if end_x_iter, break; end
         grad{1} = conv2(res{1}, rot90(k,2), 'same');
         grad{2} = conv2(res{2}, rot90(k,2), 'same');
         delta = opt_step_size(res,grad,k);
         
         ls_iter = 1; success = false;
         while ~success && ls_iter <= max_ls_iter,
            vf1 = @() x{1} - delta*grad{1};
            vf2 = @() x{2} - delta*grad{2};
            
            xTr = proj({vf1(),vf2()});
            
            % optimality test
            if ls_iter ==1 && optTol > 0,
               optCond = ...
                  norm((xTr{1}-x{1})/delta,'fro')^2 + ...
                  norm((xTr{2}-x{2})/delta,'fro')^2;
               if optCond < optTol, end_x_iter = true; end
            end
            
            resTr{1} = conv2(xTr{1}, k, 'same') - y{1};
            resTr{2} = conv2(xTr{2}, k, 'same') - y{2};
            
            lTrCost = (1/2)*(norm(resTr{1},'fro')^2 + norm(resTr{2},'fro')^2);
            pTrCost = conFun(xTr);
            
            
            decr = sum(sum(grad{1}.*(xTr{1}-x{1}) + grad{2}.*(xTr{2}-x{2})));
            sdFac = -Inf;%1e-4;
            if lTrCost - lcost(ti-1) <= sdFac*decr,
               success = true;
               lcost(ti) = lTrCost;
               pcost(ti) = pTrCost;
               x = xTr;
               res = resTr;
               ti = ti+1;
            else
               disp('backtrack');
               delta = delta/4;
            end
            
            if ls_iter > 2
               kHat = fft2(k,size(k,1)*2,size(k,2)*2);
               delta = 0.99/max(vec(abs(kHat).^2));
            end
            ls_iter = ls_iter+1;
         end
         predL1L2 = L12g({res{1}+y{1},res{2}+y{2}})^2;
         pa_ratio = predL1L2/ypcost;
         
         if ls_iter > max_ls_iter, warning('**LINE SEARCH FAILURE**'); end
         if (pTrCost / t > 1.00001), warning('** PROJECTION VIOLATION **'); end
         
      end
            
      %%%
      %%% adjust sparsity constraint (t)
      %%%
      
      if mod(iter,stageInd) == 0 && iter >= burnIn,
         if pa_ratio < pa_target,
            DFacMax = 100;
            tFac = min(DFacMax,tJumpFac*tFac);
         end
         
         t = min(numel(x{1})-2,max(4,tFac*ypcost));
         
         proj = @(x) vec22gc(vec_prox(@top_k_hard,gc2vec2(x),round(t)));
      end
      
      %%%
      %%% k update: SPG with sparse simplex projection     
      %%%
      
      size_k = size(k);
      % precompute RHS
      for iii = 1:length(x)
         flipX{iii} = rot90(x{iii},2);
      end
      kerFunObj = @(k) bdc_obj(reshape(k,size_k), x, flipX, y,'psf');
      kerFunProj = @(k) projectSimplex(k(:));
      spgOpts.maxIter = k_in_iter;
      spgOpts.verbose = 0;
      spgOpts.curvilinear = 1;
      spgOpts.optRelTol = 1e-2;
      [spg_k,spg_f,spg_funEvals] = minConf_SPG(kerFunObj,k(:),kerFunProj,spgOpts);
      
      kLog(:,:,iter-1) = k;
      k = reshape(spg_k,size_k);
      
      ektime = toc(t0k);
      edtime = toc(t0d);
      if showFigs && (ektime > timePerKernDisp),
         figure(1); colormap('gray');
         subplot(311);
         imagesc(k.^0.5); colorbar;
         title('kernel');
         %set(gca,'CLim',[0 1]);
         axis image;
         subplot(312);
         imagesc(max(0,k./kLog(:,:,max(1,iter-1))-1));
         set(gca,'CLim',[0 0.2]); colorbar;
         title('kernel frac change (+)');
         axis image;
         subplot(313);
         imagesc(max(0,-k./kLog(:,:,max(1,iter-1))+1));
         title('kernel frac change (-)');
         set(gca,'CLim',[0 0.2]); colorbar;
         axis image;
         
         drawnow;
         t0k = tic;
      end
      if showFigs && (edtime > timePerDataDisp),
         xdisp = reshape(normd([vec(x{1}), vec(x{2})],2,2),size(x{1}));
         figure(2);
         ax(1)=subplot(121);
         imagesc(xdisp);
         set(gca,'CLim',[0,prctile(xdisp(:),100)]);
         axis image;
         title('gradient norm');
         ax(2)=subplot(122);
         imagesc(...
            reshape(normd([vec(res{1}), vec(res{2})],2,2),size(res{1})));
         axis image;
         title('residual');
         linkaxes(ax);
         t0d = tic;
         clear xdisp;
         
      end
      
      fprintf('%5d %13.4e %13.4e %13.4f %13.4f\n', ...
         iter,lcost(ti-1),t,t/ypcost,pa_ratio);
      
   end;
   
   
end


function delta = opt_step_size(res,grad,k)
cgrad{1} = conv2(grad{1},k,'same');
cgrad{2} = conv2(grad{2},k,'same');
vtv = sum(sum(cgrad{1}.^2 + cgrad{2}.^2));
rtv = sum(sum(cgrad{1}.*res{1} + cgrad{2}.*res{2}));
delta = rtv/vtv;
end
