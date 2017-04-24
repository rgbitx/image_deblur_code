function [mu,sigma,weight,log_likelihood]=GaussianMixtures1D(x,nComponents)

% Author: Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology


%
% EM implementation for mixture of 2-D Gaussians.
%
% Is fixed to 2-D data I'm afriad.
%
% Inputs
% ------
% x - x-coordinate data (1 x nPoints float)
% y - y-coordinate data (1 x nPoints float)
% nComponents - # of components in mixure model (1 x 1 int)
%
% Outputs
% -------
% mu - final estimate of means (2 x nComponents float)
% sigma - final estimate of covariance matrices (2 x 2 x nComponents float)
% weight - components weights (1 x nComponents float)
% log_likelihood - progression of log-likelihood (1 x nIterations float)
%
% Example
% -------
% [x,y,gt_mu,gt_sigma,gt_weight,bg_points]=CreateMixture4(1000,3,10,1,0);
% [mu,sigma,weight,log_likelihood]=GaussianMixtures(x,y,3); 
%
%
% Rob Fergus 24/5/05 for VGG reading group
%

% Control parameters
DEBUG                       = 0;    %% Turns plotting off/on
MAX_ITERATIONS              = 100;  %% Hard upper limit on # iterations
RAND_INIT_COVARIANCE_MATRIX = 0;    %% 0 - use isotropic Gaussians, 1 -
                                    %% use randomly generated coviance matrices
LIKELIHOOD_CHANGE_THRESHOLD = 1e-5; %% termination criterion (min
                                    %% likelihood change from one iteration
                                    %% to next)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Get data so it is in column format.
if (size(x,1)>size(x,2))
   x=x';
end

nPoints = size(x,2);
cols={'r' 'g' 'b' 'c' 'm' 'y' 'k'};
ell_handle = zeros(1,nComponents);
delta_lh = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise start values for pi, mu and Sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for a=1:nComponents
   
   %%% Get mean of data and peturb a bit
   xinit = mean(x) + (rand-0.5);
   
   %%% Set 
   mu(1,a)=0; %xinit;

   if RAND_INIT_COVARIANCE_MATRIX
      %%% Generate random matrix
      sigma(:,:,a)=eye(2)*0.5 + RandCovarianceMatrix(2,1);   
   else
      %%% Just set to identitiy matrix
      sigma(:,:,a)=(1e6-(rand(1)*1e6))*eye(1); 
   end
   
   %%% Equal weight to each component initially
   weight(a)=1/nComponents;  
   
end

sigma(:,:,1) = 1e6;

%%%% Create empty matrix for responsibilities
%
resp        = zeros(nComponents,nPoints);
likelihoods = zeros(nComponents,nPoints);	
iteration   = 1;

for iteration = 1:MAX_ITERATIONS
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % E-Step
   %
   % Calc. new responsibilities
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for c=1:nComponents
  
      %%% Compute normaliser for Gaussian
      normaliser = 1/((2 * pi)^0.5) * det(sigma(:,:,c))^-0.5;
      tmp = mu(:,c) * ones(1,nPoints);
   
      %%% Precompute (x-mu) bit (since it is used twice)
      offset_data = [x] - tmp;
     
      %%% Compute exponent part of Gaussian
      exponent = sum(offset_data .*  ( inv(sigma(:,:,c)) * offset_data ),1);
      
      %%% Now compute responsibilites (unnormalised)
      likelihoods(c,:) =  weight(c) * normaliser * exp( -0.5 * exponent );
      
   end

   %%% Compute log-likelihood.....
   %%% Log-likelihood 1/(# points) * \sum_points log \sum_components
   %%% p(point x | component c) 
   log_likelihood(iteration) = mean(log(sum(likelihoods,1)));
   
   if (iteration>1)
      delta_lh(iteration) = log_likelihood(iteration)-log_likelihood(iteration-1);
      fprintf('Iteration: %d, log-likelihood: %f delta-ll: %f\n',iteration,log_likelihood(iteration),delta_lh(iteration));
   else
     fprintf('Iteration: %d, log-likelihood: %f\n',iteration,log_likelihood(iteration));
   end
   
   %%% Now normalise likelihoods to get repsonsibiltiies: p^{t+1}(z|x,theta^t)
   resp_total = sum(likelihoods,1);
   
   for c=1:nComponents % and normalise over the mixture
   	resp(c,:) = likelihoods(c,:) ./ resp_total;
   end      
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % M-step
   %
   % Update parameters
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   for c=1:nComponents
   
      %%% Get total responsibility for this component
      %%% \sum_points p(z=c|x,\theta^t)
      total_resp_this_component = sum(resp(c,:));
      
      % Update component weights
      weight(c) = total_resp_this_component / nPoints;
      
      
      % Calc. new mean
      mu(1,c) = 0;
            
      % and new covariance matrices
      % use trick to avoid loop over all points...
      u=[sqrt(resp(c,:))];   
      
      %%% note we use updated versions of mu, not old ones....
      new_offset_data = [x] - ( mu(:,c) * ones(1,nPoints) );
      
      %%% now actually compute sigma. note lack of look, just clever
      %%% inner products
      sigma(:,:,c) = (( u.* new_offset_data ) * ( u .* new_offset_data)') / total_resp_this_component;
      
      sigma(:,:,c) = sigma(:,:,c) + 1e-5;   
   end

   if (iteration<10)
     sigma(:,:,1) = 1e6;
   end
   
   %%%%% Termination section
   if (delta_lh(iteration)<LIKELIHOOD_CHANGE_THRESHOLD)
      break;
   end
      
   if DEBUG
      
      figure(1); clf;
      plot(x,y,'b.');
      title(['Iteration: ',num2str(iteration)]);
      for c=1:nComponents
	    lw = round( 10 * weight(c) ) + 1;
	    [ell_handle(c),xy,L,l,th] = draw_ellipse(mu(:,c),sigma(:,:,c),cols{rem(c-1,7)+1},lw,[]);
      end
      
      pause
   end
      
end


if DEBUG
   figure(2); clf;
   plot([1:iteration],(log_likelihood),'b.-');
   xlabel('Iteration');
   ylabel('Log-likelihood');
   title('Evolution of log-likelihood');
end
