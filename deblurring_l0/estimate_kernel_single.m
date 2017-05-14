function [K,L,opts] = estimate_kernel_single(B,K,opts)
% derivative filters
B = double(B);
if max(B(:)) > 2
    B = B/255;
end
for iter = 1:opts.xk_iter      
   %% get latent image    
   
   L = L0Restoration_mul(B,K,opts.lambda); 
   
   %% get kernel   
   
   [K,opts] = estimate_psf_gradient_mul(B,L,opts,K);  
   
   fprintf('%d iterations is done\n', iter);
   
    if opts.verbose
        figure(1000),imshow(mat2gray(K(:,:)));
        figure(2000),imshow(mat2gray(L));
        title(sprintf('iter=%d',iter));
        pause(0.01);    
    end
    % update lambda
    if mod(iter,opts.iter_step) == 0
        opts.lambda = max(opts.min_lambda,opts.lambda/1.1);        
    end%if
end%iter

% figure, imshow(K);
% title('kernel')

end