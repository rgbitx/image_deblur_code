function [ensemble,grad]=train_ensemble_evidence6(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12)

% Author: James Miskin, adapted by Rob Fergus
% Version: 1.0, distribution code.
% Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
% Copyright 2006, Massachusetts Institute of Technology

%
% Code from James Miskin, adapted by R.Fergus for use in SIGGRAPH 06
% paper 'Removing Camera Shake from a Single Photograph' by R. Fergus,
% B. Singh, A. Hertzmann, S. T. Roweis and W. T. Freeman
%

% Original code available at:
% http://www.inference.phy.cam.ac.uk/jwm1003/train_ensemble.tar.gz
% For details of code operation, please see the following publications:
%
% @incollection{Miskin2000,
% author = {James Miskin and David J. C. MacKay},
% title = {Ensemble {L}earning for {B}lind {I}mage {S}eparation and {D}econvolution},
% booktitle = {Adv.~in {I}ndependent {C}omponent {A}nalysis},
% editor = {M. Girolani},
% publisher = {Springer-Verlag},
% year = {2000}
% }
%
% @phdthesis {Miskin2000a,
%        author = "Miskin, James W.",
%        title = {Ensemble Learning for Independent Component Analysis},
%        school = {University of Cambridge},
%        year = {2000},
%        month = {December},
%        url = {http://www.inference.phy.cam.ac.uk/jwm1003/} }
  
%=========================================================================
% Evaluates evidence for the current ensemble.
%   step_len       - Size of step to take
%   dimensions     - Class description matrix
%   opt_func       - Function that returns the optimum parameters for
%                    distribution5
%   ensemble       - Current ensemble
%   direction      - Direction to search along in ensemble space
%   state          - State variable 
%                      1 - Don't update priors, noise
%                      2 - Don't update noise
%                      3 - Update all
%   P1-P10         - Variables that are passed through to opt_func
%   P9 - prior structure with fields pi & gamma
%   P10 - fft mode. 1 = use fft convolutions, 0 = use conv2
%   P11 - blur mask (vector size of blur) giving weighting on the blur
%   elements
%   P12 - image mask (vector size of image) giving weighting on variance
%   of image elements - used to alter saturation  
%=========================================================================

%Set string to evaluate optimum
%optstr = [opt_func];
%optstr=[optstr, '(dimensions,ensemble'];
%for i=1:nargin - 6
%  optstr = [optstr,',P',int2str(i)];
%end
%optstr = [optstr, ')'];
%train_blind_deconv2(dimensions,ensemble,P1,P2,P3,P4,P5,P6,P7,P8)

%Step to the test point
ensemble.x1    =ensemble.x1    +step_len*direction.x1;
ensemble.x2    =abs(ensemble.x2    +step_len*direction.x2);
ensemble.b_x_2 =abs(ensemble.b_x_2 +step_len*direction.b_x_2);
ensemble.ba_x_2=abs(ensemble.ba_x_2+step_len*direction.ba_x_2);
ensemble.pi_x_2=abs(ensemble.pi_x_2+step_len*direction.pi_x_2);

%Make a copy of the direction
grad=direction;

%Check for valid ensemble
if (any(ensemble.x2<0) | any(ensemble.b_x_2<0) | any(ensemble.ba_x_2<0) | any(ensemble.pi_x_2<0))
  ensemble.D_val=Inf;
  grad.x1(:)=NaN;
  grad.x2(:)=NaN;
  grad.b_x_2(:)=NaN;
  grad.ba_x_2(:)=NaN;
  grad.pi_x_2(:)=NaN;
  return;
end

%Find evidence by summing over the components
ptr=0;
for c=1:size(dimensions,1)

    %Update the expectations based on the updated Q(x)
    cx1 =train_ensemble_get(c,dimensions,ensemble.x1);
    cx2 =train_ensemble_get(c,dimensions,ensemble.x2);
   
    [cHx,cmx,cmx2]=train_ensemble_rectified5(cx1,cx2,dimensions(c,4)); 

  %Find optimum weights
  c_pi_x_2=reshape(ensemble.pi_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)),dimensions(c,1),dimensions(c,3));
  c_b_x_2=reshape(ensemble.b_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)),dimensions(c,1),dimensions(c,3));
  c_ba_x_2=reshape(ensemble.ba_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)),dimensions(c,1),dimensions(c,3));
  c_log_lambda_x=zeros(dimensions(c,1:3));
  
  if (dimensions(c,3)>1)
    if (dimensions(c,4)==0 | dimensions(c,4)==2)
      %Gaussian or Rectified Gaussian prior
      for alpha=1:dimensions(c,3)
	for k=1:dimensions(c,1)
	  c_log_lambda_x(k,:,alpha) = (log(c_pi_x_2(k,alpha))-0.5/ c_pi_x_2(k,alpha)+0.5*log(c_ba_x_2(k,alpha))-0.25/c_b_x_2(k,alpha)-0.5*cmx2(k,:)*c_ba_x_2(k,alpha));
        end
      end
    elseif (dimensions(c,4)==1)
      %Exponential prior
      for alpha=1:dimensions(c,3)
	if isempty('P11')
          keyboard
          for k=1:dimensions(c,1)
            c_log_lambda_x(k,:,alpha) = (log(c_pi_x_2(k,alpha))-0.5/c_pi_x_2(k,alpha)+log(c_ba_x_2(k,alpha))-0.5/c_b_x_2(k,alpha)-cmx(k,:)*c_ba_x_2(k,alpha)) + P11(alpha,:);
          end
        else
          for k=1:dimensions(c,1)
            c_log_lambda_x(k,:,alpha) = (log(c_pi_x_2(k,alpha))-0.5/c_pi_x_2(k,alpha)+log(c_ba_x_2(k,alpha))-0.5/c_b_x_2(k,alpha)-cmx(k,:)*c_ba_x_2(k,alpha));
          end
        end
      end
    elseif (dimensions(c,4)==3)
      %Discrete prior      
    elseif  (dimensions(c,4)==4)
      %Exponential prior
      for alpha=1:dimensions(c,3)
	for k=1:dimensions(c,1)
	  c_log_lambda_x(k,:,alpha)=log(c_pi_x_2(k,alpha))-0.5/c_pi_x_2(k,alpha)+0.5*log(c_ba_x_2(k,alpha))-0.25/c_b_x_2(k,alpha)-0.5*abs(cmx(k,:))*c_ba_x_2(k,alpha);
	end
      end
    end
	
    %Now normalise c_log_lambda_x
    max_c_log_lambda_x=max(c_log_lambda_x,[],3);
    for alpha=1:dimensions(c,3)
      c_log_lambda_x(:,:,alpha)=c_log_lambda_x(:,:,alpha)-max_c_log_lambda_x;
    end
    log_sum_lambda_x=log(sum(exp(c_log_lambda_x),3));
    for alpha=1:dimensions(c,3)
      c_log_lambda_x(:,:,alpha)=c_log_lambda_x(:,:,alpha)-log_sum_lambda_x;
    end 
  end
  
  %Find the optimum parameters for the prior parameters
  if (dimensions(c,4)==0 | dimensions(c,4)==2)
    %(Rectified) Gaussian prior for P(x)
    opt_c_b_x_2=ensemble.b_x+reshape(sum(exp(c_log_lambda_x),2),dimensions(c,1),dimensions(c,3))/2;
    opt_c_pi_x_2=ensemble.pi_x+reshape(sum(exp(c_log_lambda_x),2),dimensions(c,1),dimensions(c,3));
    opt_c_ba_x_2=zeros(dimensions(c,1),dimensions(c,3));
    for alpha=1:dimensions(c,3)
      opt_c_ba_x_2(:,alpha)=opt_c_b_x_2(:,alpha)./(ensemble.a_x+sum(exp(c_log_lambda_x(:,:,alpha)).*cmx2,2)/2);
    end
  elseif (dimensions(c,4)==1)
    %Exponential prior for P(x)
    opt_c_b_x_2=ensemble.b_x+reshape(sum(exp(c_log_lambda_x),2),dimensions(c,1),dimensions(c,3));
    opt_c_pi_x_2=ensemble.pi_x+reshape(sum(exp(c_log_lambda_x),2),dimensions(c,1),dimensions(c,3));
    opt_c_ba_x_2=zeros(dimensions(c,1),dimensions(c,3));
    for alpha=1:dimensions(c,3)
      opt_c_ba_x_2(:,alpha)=opt_c_b_x_2(:,alpha)./(ensemble.a_x+sum(exp(c_log_lambda_x(:,:,alpha)).*cmx,2));
    end
  elseif ((dimensions(c,4)==3) | (dimensions(c,4)==4))
    %Discrete prior (not learning any parameters of prior)
    opt_c_b_x_2=ensemble.b_x*ones(dimensions(c,1),dimensions(c,3));
    opt_c_ba_x_2=ensemble.b_x/ensemble.a_x*ones(dimensions(c,1),dimensions(c,3));
    opt_c_pi_x_2=ensemble.pi_x*ones(dimensions(c,1),dimensions(c,3));
  elseif (dimensions(c,4)==5)
    %Laplacian prior for P(x)
    opt_c_b_x_2=ensemble.b_x+reshape(sum(exp(c_log_lambda_x),2),dimensions(c,1),dimensions(c,3));
    opt_c_pi_x_2=ensemble.pi_x+reshape(sum(exp(c_log_lambda_x),2),dimensions(c,1),dimensions(c,3));
    opt_c_ba_x_2=zeros(dimensions(c,1),dimensions(c,3));
    for alpha=1:dimensions(c,3)
      opt_c_ba_x_2(:,alpha)=opt_c_b_x_2(:,alpha)./(ensemble.a_x+sum(exp(c_log_lambda_x(:,:,alpha)).*cmx,2));
    end
  end

  %%%% Manual over-ride for priors
  if (dimensions(c,5)>0)
 
    if (c==3) %% Image prior
      opt_c_pi_x_2 = P9.pi(1:dimensions(c,1),:) * 1e3;
      opt_c_ba_x_2 = P9.gamma(1:dimensions(c,1),:);
      opt_c_b_x_2 = ones(dimensions(c,1),dimensions(c,3)) * 1e-3; 
    elseif (c==2) %% Blur prior
  
      opt_c_ba_x_2 = [ 5.1143e3 5.0064e+03 173.8885 50.6538 ];
      opt_c_b_x_2 = [787.8988 201.7349 236.1948 143.1756];
    else
      error foo
    end
      
  else
    % Leave alone
  end   

  
  %Optimum Q(x)
  opt_cx1=zeros(size(cx1));
  opt_cx2=zeros(size(cx2));
  if (dimensions(c,4)==0 | dimensions(c,4)==2)
    for alpha=1:dimensions(c,3)
      for k=1:dimensions(c,1)
	opt_cx2(k,:)=opt_cx2(k,:)+c_ba_x_2(k,alpha)*exp(c_log_lambda_x(k,:,alpha));
      end
    end
  elseif ((dimensions(c,4)==1)| dimensions(c,4)==4)
    for alpha=1:dimensions(c,3)
      for k=1:dimensions(c,1)
	opt_cx1(k,:)=opt_cx1(k,:)-c_ba_x_2(k,alpha)*exp(c_log_lambda_x(k,:,alpha));
      end
    end
  end

  %Evaluate the KL divergence for the distributions trained in this function
  ensemble.D_x(sum(dimensions(1:c-1,1))+1:sum(dimensions(1:c,1)))=sum(cHx,2)+...
      sum(gammaln(ensemble.b_x)-ensemble.b_x*log(ensemble.a_x)-gammaln(c_b_x_2)+c_b_x_2.*log(c_b_x_2./c_ba_x_2)+(c_b_x_2-opt_c_b_x_2).*(log(c_ba_x_2)-0.5./c_b_x_2)+(opt_c_b_x_2./opt_c_ba_x_2-c_b_x_2./c_ba_x_2).*c_ba_x_2,2)+...
      sum(gammaln(ensemble.pi_x)-gammaln(c_pi_x_2)+(c_pi_x_2-opt_c_pi_x_2).*(log(c_pi_x_2)-0.5./c_pi_x_2),2)+...
      (-gammaln(dimensions(c,3)*ensemble.pi_x)+gammaln(sum(c_pi_x_2,2))+sum(c_pi_x_2-opt_c_pi_x_2,2).*(-log(sum(c_pi_x_2,2))+0.5./sum(c_pi_x_2,2)))+...
      sum(sum(c_log_lambda_x.*exp(c_log_lambda_x),2),3);

  %Store updated parameters and optimum distributions
  ensemble.log_lambda_x=train_ensemble_put_lambda(c,dimensions,ensemble.log_lambda_x,c_log_lambda_x);
  ensemble.mx=train_ensemble_put(c,dimensions,ensemble.mx,cmx);
  ensemble.mx2=train_ensemble_put(c,dimensions,ensemble.mx2,cmx2);
  ensemble.opt_ba_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)) = reshape(opt_c_ba_x_2,1,dimensions(c,1)*dimensions(c,3));
  ensemble.opt_b_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3))  = reshape(opt_c_b_x_2 ,1,dimensions(c,1)*dimensions(c,3));
  ensemble.opt_pi_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)) = reshape(opt_c_pi_x_2,1,dimensions(c,1)*dimensions(c,3));  
  
  grad.x1=train_ensemble_put(c,dimensions,grad.x1,opt_cx1);
  grad.x2=train_ensemble_put(c,dimensions,grad.x2,opt_cx2);
  
  %Increment the pointer
  ptr=ptr+dimensions(c,1)*dimensions(c,3);
 
end

%Only set the gradient for the priors if in state 2 or 3
if (state>=2)
  grad.pi_x_2=ensemble.opt_pi_x_2;
  grad.b_x_2 =ensemble.opt_b_x_2;
  grad.ba_x_2=ensemble.opt_ba_x_2;
else
  grad.pi_x_2=ensemble.pi_x_2;
  grad.b_x_2 =ensemble.b_x_2;
  grad.ba_x_2=ensemble.ba_x_2;
end

%Q(gamma)
if (P10)
  %%% fft mode
  [dx1,dx2,rerror,data_points]=train_blind_deconv(dimensions,ensemble,P1,P2,P3,P4,P5,P6,P7,P8);
else
  [dx1,dx2,rerror,data_points]=train_blind_deconv2(dimensions,ensemble,P1,P2,P3,P4,P5,P6,P7,P8);
end
ensemble.b_sigma_2=ensemble.b_sigma+data_points/2;
ensemble.opt_ba_sigma_2=ensemble.b_sigma_2/(ensemble.a_sigma+rerror/2);
%Set to the optimum if in the final state
if (state==3)
  ensemble.ba_sigma_2=ensemble.opt_ba_sigma_2;
end

%Q(x)
grad.x1=grad.x1+ensemble.ba_sigma_2*dx1;
grad.x2=grad.x2+ensemble.ba_sigma_2*dx2;

%Evaluate the KL divergence
ensemble.D_x(sum(dimensions(:,1))+1)=gammaln(ensemble.b_sigma)-gammaln(ensemble.b_sigma_2)-ensemble.b_sigma*log(ensemble.a_sigma)+ensemble.b_sigma_2*log(ensemble.b_sigma_2/ensemble.ba_sigma_2)+(ensemble.b_sigma_2/ensemble.opt_ba_sigma_2-ensemble.b_sigma_2/ensemble.ba_sigma_2)*ensemble.ba_sigma_2+data_points*log(2*pi)/2;

%Normalise from log_e to bits per data point
ensemble.D_x=ensemble.D_x/data_points*log2(exp(1));
if (isnan(ensemble.D_val))
  ensemble.D_val=inf;
end

%Evaluate the KL divergence
ensemble.D_val=sum(ensemble.D_x);

%Find the direction
grad.x1=grad.x1-ensemble.x1;
grad.x2=grad.x2-ensemble.x2;
grad.b_x_2=grad.b_x_2-ensemble.b_x_2;
grad.ba_x_2=grad.ba_x_2-ensemble.ba_x_2;
grad.pi_x_2=grad.pi_x_2-ensemble.pi_x_2;

%Don't train classes if options spectify clamping
ptr=0;
for c=1:size(dimensions,1)
  if (~dimensions(c,6))
    %%% Comment out line below & set dimensions(c,6) to 0 for MAP
    %%% estimate 
    grad.x1(ptr+1:ptr+dimensions(c,1)*dimensions(c,2))=0;
    grad.x2(ptr+1:ptr+dimensions(c,1)*dimensions(c,2))=0;
  end
  ptr=ptr+dimensions(c,1)*dimensions(c,2);
end

%% Image masking to correct for saturation
ptr=0;
for c=1:size(dimensions,1)
  if (c==3) %%% image only
    %grad.x1(ptr+1:ptr+dimensions(c,1)*dimensions(c,2))=grad.x1(ptr+1:ptr+dimensions(c,1)*dimensions(c,2)) .* P12';

%    grad.x2(ptr+1:ptr+dimensions(c,1)*dimensions(c,2))=grad.x2(ptr+1:ptr+dimensions(c,1)*dimensions(c,2)) .* P12';
  end   
  ptr=ptr+dimensions(c,1)*dimensions(c,2);
end
