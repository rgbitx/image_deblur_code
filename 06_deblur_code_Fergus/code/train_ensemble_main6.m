function [ensemble,D_log,gamma_log]=train_ensemble_main6(dimensions,initial_x1,initial_x2,opt_func,text,options,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12)

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
% Performs ensemble learning to model a data set.
%   dimensions     - Class description matrix
%   initial_x1,x2  - Initial parameters of the ensemble, can be []
%   opt_func       - Function that returns the optimum parameters for
%                    distributions
%   text           - Text that is displayed at each iteration
%   options        - Options to control the training
%   P1-P10         - Variables that are passed through to opt_func
%=========================================================================
% Performs ensemble fitting of posterior over parameters by a separable
% mixture distribution. The problem is described by the dimensions array
% which stores a row for each type of data in the problem. Each row has
% the following columns (for row i)
%  1 - Number of rows of x_i (each row has a different prior)
%  2 - Number of columns of x_i
%  3 - Number of components to use for the mixture prior
%  4 - Type of prior to use for x_i
%        0 - Gaussian
%        1 - Exponential
%        2 - Rectified Gaussian
%
% When training, opt_func is called to evaluate the optimal
% distributions. It is called as
%   opt_func(dimensions,ensemble,P1,P2,P3,...)
% The routine must return the number of data points, the error in the
% reconstruction and the optimum values for the parameters of the
% distributions (without including the prior). See the example models to
% see how this works. The priors that have been implemented so far result
% in a posterior is of the form
%    log Q(x) = -x2*x/2+x1*x + constant
% so the model file returns the optimum x1,x2 parameters given the
% data. The library routines add in the terms from the priors.
%
% options controls the training
%   options(1) - Convergence criteria, the algorithm is assumed to be
%                converging if the improvement at each iteration is less
%                than the criteria.
%   options(2) - non-zero if the cost function along the search direction
%                is to be displayed at each iteration (slow).
%   options(3) - Initial noise variance to be assumed.
%   options(4) - Non-zero if priors are to be reinitialised whenever the
%                noise changes
%   options(5) - Non-zero if any components that are switched off (ie
%                tend to prior) are to be re-randomised when noise
%                changes, this may be an advantage if the noise is
%                initially assumed to be large.
%=========================================================================

%Create array to store training logs
Niter=options(6);
D_log=NaN*ones(2,Niter);
gamma_log=NaN*ones(1,Niter);

if ~exist('P11')
  P11 = [];
end

if ~exist('P12')
  P12 = [];
end

%Create the ensemble parameters
ensemble=struct('x1' ,1e4*randn(1,dimensions(:,1)'*dimensions(:,2)).*ceil(rand(1,dimensions(:,1)'*dimensions(:,2))*2),...
		'x2',1e4*ones(1,dimensions(:,1)'*dimensions(:,2)),...
		'mx' ,zeros(1,dimensions(:,1)'*dimensions(:,2)),...
		'mx2',zeros(1,dimensions(:,1)'*dimensions(:,2)),...
		'log_lambda_x',zeros(1,sum(dimensions(:,1).*dimensions(:,2).*dimensions(:,3))),...
		'pi_x',1,...
		'pi_x_2',ones(1,dimensions(:,1)'*dimensions(:,3)),...
		'opt_pi_x_2',ones(1,dimensions(:,1)'*dimensions(:,3)),...
		'a_x',1e-3,...
		'ba_x_2' ,ones(1,dimensions(:,1)'*dimensions(:,3)),...
		'opt_ba_x_2' ,ones(1,dimensions(:,1)'*dimensions(:,3)),...
		'b_x',1e-3,...
		'b_x_2' ,ones(1,dimensions(:,1)'*dimensions(:,3)),...
		'opt_b_x_2' ,ones(1,dimensions(:,1)'*dimensions(:,3)),...
		'a_sigma',1e-3,...
		'ba_sigma_2',0,...
		'opt_ba_sigma_2',0,...
		'b_sigma',1e-3,...
		'b_sigma_2',0,...
		'D_val',0,...
		'D_x',zeros(sum(dimensions(:,1))+1,1));

direction=struct('x1' ,zeros(1,dimensions(:,1)'*dimensions(:,2)),...
		 'x2' ,zeros(1,dimensions(:,1)'*dimensions(:,2)),...
		 'pi_x_2',zeros(1,dimensions(:,1)'*dimensions(:,3)),...
		 'ba_x_2' ,zeros(1,dimensions(:,1)'*dimensions(:,3)),...
		 'b_x_2' ,zeros(1,dimensions(:,1)'*dimensions(:,3)));
	 
%Set the initial values of the distributions
if (~isempty(initial_x1))
  ensemble.x1    =initial_x1;
end
if (~isempty(initial_x2))
  ensemble.x2    =initial_x2;
end



%Make sure no parameters are set to zero initially
ensemble.x1=ensemble.x1+1e-16*(ensemble.x1==0);
ensemble.x2=ensemble.x2+1e-16*(ensemble.x2==0);

%Set string to evaluate evidence
%evidencestr='train_ensemble_evidence5(step_len,dimensions,opt_func,ensemble,direction,state';
%for i=1:nargin - 6
%  evidencestr = [evidencestr,',P',int2str(i)];
%end
%evidencestr = [evidencestr, ')'];
%train_ensemble_evidence5(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8)
%keyboard

%Details of the state
state=1;
last_change_iter=0;
alpha=1.0;
beta=0.9;

%Copy options
converge_criteria=options(1);
plot_step=options(2);
ensemble.ba_sigma_2=options(3).^-2;
restart_priors=options(4);
restart_switched_off=options(5);

%Iterate algorithm to improve the approximating ensemble
oD_val=NaN;
for iter=1:Niter

  %Re-evaluate after a state change
  if (iter==last_change_iter+1)
    if (state<3)
      %Evaluate the model before updating the priors
      step_len=0.0;
      
      [ensemble,grad]=train_ensemble_evidence6(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12);
      direction=grad;

      %Set the priors to their optimal values
      ensemble.pi_x_2=ensemble.opt_pi_x_2;
      ensemble.b_x_2 =ensemble.opt_b_x_2;
      ensemble.ba_x_2=ensemble.opt_ba_x_2;
            
      %Re-initialise priors that have coallesced
      ptr=0;
      for c=1:size(dimensions,1)
	%Only need to check if there is a mixture
	if (dimensions(c,3)>1)
	  c_pi_x_2=reshape(ensemble.pi_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)),dimensions(c,1),dimensions(c,3));
	  c_b_x_2=reshape(ensemble.b_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)),dimensions(c,1),dimensions(c,3));
	  c_ba_x_2=reshape(ensemble.ba_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)),dimensions(c,1),dimensions(c,3));
      
	  for k=1:dimensions(c,1)
	    %Check for each distribution type for coallescence
	    if (dimensions(c,4)==0 | dimensions(c,4)==1 | dimensions(c,4)==2)
	      %Gaussian or Exponential
	      sorted_scales=sort(c_ba_x_2(k,:));
	      if (restart_priors | any(c_pi_x_2(k,:)<ensemble.pi_x+1/dimensions(c,3)) | any(sorted_scales(2:end)<1.5*sorted_scales(1:end-1)))
		%If any component has approximately zero weight or if any
                %two have approximately the same scale, they should be
                %restarted
		mean_scale=sum(c_b_x_2(k,:)./c_ba_x_2(k,:))/sum(c_b_x_2(k,:));
		c_pi_x_2(k,:)=ensemble.pi_x+dimensions(c,2)/dimensions(c,3);
		c_b_x_2(k,:) =ensemble.b_x +dimensions(c,2)/dimensions(c,3);
		c_ba_x_2(k,:)=c_b_x_2(k,:)./(ensemble.a_x+0.5*(1:dimensions(c,3))*mean_scale*dimensions(c,2)/dimensions(c,3));
              end		
	    elseif (dimensions(c,4)==3 | dimensions(c,4)==4)
	      %Discrete so no prior properties
	    end
	  end
	   
	  %Store the new parameters
	  ensemble.pi_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3))=reshape(c_pi_x_2,1,dimensions(c,1)*dimensions(c,3));
	  ensemble.b_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3)) =reshape(c_b_x_2,1,dimensions(c,1)*dimensions(c,3));
	  ensemble.ba_x_2(ptr+1:ptr+dimensions(c,1)*dimensions(c,3))=reshape(c_ba_x_2,1,dimensions(c,1)*dimensions(c,3));
	end
	
	%Increment the pointer
	ptr=ptr+dimensions(c,1)*dimensions(c,3);
      end
    end
    
    %Re-evaluate the evidence and the search direction
    step_len=0.0;
         
    [ensemble,grad]=train_ensemble_evidence6(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12);
      
    direction.x1=alpha*grad.x1;
    direction.x2=alpha*grad.x2;
    direction.b_x_2=alpha*grad.b_x_2;
    direction.ba_x_2=alpha*grad.ba_x_2;
    direction.pi_x_2=alpha*grad.pi_x_2;
    step_len=1.0;
  end
  
  if (plot_step)
    ostep_len=step_len;
    x_vals=(0:0.1:2)*step_len;
    f_vals=zeros(size(x_vals));
    for i=1:prod(size(x_vals))
      step_len=x_vals(i);
      [pensemble,pgrad]=train_ensemble_evidence6(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12);
      f_vals(i)=pensemble.D_val;
    end
  %  figure(10)
  %  plot(x_vals,f_vals)
  %  v=axis;
  %  line([ostep_len ostep_len],v(3:4))
  %  step_len=ostep_len;
  %  pause
  end

 % drawnow;

  [tensemble,tgrad]=train_ensemble_evidence6(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12);
  success=(tensemble.D_val<ensemble.D_val);
  if (tensemble.D_val>ensemble.D_val)
    direction.x1=alpha*grad.x1;
    direction.x2=alpha*grad.x2;
    direction.b_x_2=alpha*grad.b_x_2;
    direction.ba_x_2=alpha*grad.ba_x_2;
    direction.pi_x_2=alpha*grad.pi_x_2;
    step_len=2*step_len;
    while (tensemble.D_val>ensemble.D_val+converge_criteria/1e4)
        %drawnow;
      step_len=0.5*step_len;
      [tensemble,tgrad]=train_ensemble_evidence6(step_len,dimensions,opt_func,ensemble,direction,state,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12);
    end
  end
  
  direction.x1=alpha*tgrad.x1+beta*direction.x1;
  direction.x2=alpha*tgrad.x2+beta*direction.x2;
  direction.b_x_2=alpha*tgrad.b_x_2+beta*direction.b_x_2;
  direction.ba_x_2=alpha*tgrad.ba_x_2+beta*direction.ba_x_2;
  direction.pi_x_2=alpha*tgrad.pi_x_2+beta*direction.pi_x_2;
  step_len=min(1,1.1*step_len);
  
  %Accept the interpolated point
  dD_val=tensemble.D_val-oD_val;
  ensemble=tensemble;
  grad=tgrad;
  
  %Store training logs
  D_log(1,iter)=ensemble.D_val;
  gamma_log(1,iter)=ensemble.ba_sigma_2;

  %Store current value
  oD_val=ensemble.D_val;
  fprintf(1, [text '  Iteration %4d  Noise=%11.6e\n'],iter,gamma_log(1,iter)^-0.5);

 
  
  %Check if algorithm has converged
  converged=0;
  if (iter>3)
    last_dD_log=D_log(1,iter-2:iter)-D_log(1,iter-3:iter-1);
    if (all(last_dD_log<converge_criteria/1e4 & last_dD_log>-converge_criteria))
      converged=1;
    end
  end
  if (converged)
    if (state==3)
      %Have finished so exit
      break;
    elseif (state==2 & ensemble.opt_ba_sigma_2<1.1*ensemble.ba_sigma_2)
      %Have converged so update noise and continue
      state=3;
    else
      %Swap between 1 and 2
      state=3-state;
    end
    last_change_iter=iter;
    
    %Run through each class and switch back on any switched off components
    if (state==1)
      %Set the noise to the optimum
      ensemble.ba_sigma_2=ensemble.opt_ba_sigma_2;
      %Randomise values for any components that have been switched off
      %(they can always switch back off again later)
      for c=1:size(dimensions,1)
	cmx=train_ensemble_get(c,dimensions,ensemble.mx);
	cmx2=train_ensemble_get(c,dimensions,ensemble.mx2);
 	if (restart_switched_off & dimensions(c,5) & dimensions(c,4)<3)
	  %Find variance scale for this class, this is the ratio of <x>^2/<x^2>
	  %For non-rectified distributions, this ratio tends to zero as the
	  %component is switched off and the posterior tends to the prior.
	  %For rectified distributions, this ratio tends to 0.6366 as the
	  %distributions tend to rectified gaussians.
	  %Therefore randomise the components if the ratio goes below 0.7
	  cx1=train_ensemble_get(c,dimensions,ensemble.x1);
	  cx2=train_ensemble_get(c,dimensions,ensemble.x2);
	  scales=mean(cmx.^2,2)./mean(cmx2,2);
	  
	  for k=1:dimensions(c,1)
	    if (scales(k)<0.7)
	      %Scale of this component is really low or C_KL is low, so
	      %randomise
	      disp(['  Reinitialising class=' int2str(c) ' k=' int2str(k)])
	      
	      if (dimensions(c,4)==0)
		%Not rectified
		cx1(k,:)=1e4*randn(1,dimensions(c,2)).*ceil(rand(1,dimensions(c,2))*2);
	      else
		%Rectified
		cx1(k,:)=1e4*abs(randn(1,dimensions(c,2)));
	      end
	      cx2(k,:)=1e4;
	    end
	  end
	  ensemble.x1=train_ensemble_put(c,dimensions,ensemble.x1,cx1);
	  ensemble.x2=train_ensemble_put(c,dimensions,ensemble.x2,cx2);
        end
      end
    end
  end
end

%Shrink the logs to use the smallest possible size
D_log=D_log(:,1:iter);
gamma_log=gamma_log(1,1:iter);











