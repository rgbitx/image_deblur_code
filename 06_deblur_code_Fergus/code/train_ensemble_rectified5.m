function [Hx,mx,mx2]=train_ensemble_rectified5(x1,x2,type)

% Author: James Miskin, adpated by Rob Fergus
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
% Evaluates expectations under the ensemble distributions. Returns the 
%   mx      = <x>
%   mx2     = <x^2>
%   Hx      = <log Q(x)> - any constants from P(x)
% The inputs are
%   x1,x2      - Parameters of the distribution
%   type       - Type of distribution
%     1 - Gaussian (so posterior is Gaussian)
%     2 - Laplacian (so posterior is Rectified Gaussian)
%     3 - Rectified Gaussian (so posterior is Rectified Gaussian)
%=========================================================================

if (type==0)
  %Gaussian
  mx=x1./x2;
  mx2=x1.^2.*x2.^-2+1./x2;
  Hx=-0.5+0.5*log(x2);
elseif (type==1)
  %Laplacian
  t=-x1./sqrt(2*x2);
  erf_table=erfcx(t);
  
  mx=(x1./x2+sqrt(2/pi./x2)./erf_table).*(t<=25)+(-1./x1+2*x2.*x1.^-3-10*x2.^2.*x1.^-5).*(t>25);
  mx2=(x1.^2.*x2.^-2+1./x2+2*x1./x2./sqrt(2*pi*x2)./erf_table).*(t<=25)+(2.*x1.^-2-10*x2.*x1.^-4+74*x2.^2.*x1.^-6).*(t>25);
  Hx=(-log(erfc(min(t,25)))+0.5*log(2*x2/pi)-0.5+x1./sqrt(2*pi*x2)./erf_table).*(t<25)+(log(-x1)-1+2*x2.*x1.^-2-15*x2.^2.*x1.^-4/2+148*x2.^3.*x1.^-6/3).*(t>=25);
elseif (type==2)
  %Rectified Gaussian
  t=-x1./sqrt(2*x2);
  erf_table=erfcx(t);
  
  mx=(x1./x2+sqrt(2/pi./x2)./erf_table).*(t<=25)+(-1./x1+2*x2.*x1.^-3-10*x2.^2.*x1.^-5).*(t>25);
  mx2=(x1.^2.*x2.^-2+1./x2+2*x1./x2./sqrt(2*pi*x2)./erf_table).*(t<=25)+(2.*x1.^-2-10*x2.*x1.^-4+74*x2.^2.*x1.^-6).*(t>25);
  Hx=(-log(erfc(min(t,25)))+0.5*log(2*x2/pi)-0.5+x1./sqrt(2*pi*x2)./erf_table).*(t<25)+(log(-x1)-1+2*x2.*x1.^-2-15*x2.^2.*x1.^-4/2+148*x2.^3.*x1.^-6/3).*(t>=25)+0.5*log(pi/2);
elseif (type==3)
  %Discrete (1,-1)
  mx=tanh(x1);
  mx2=ones(size(x1));
  Hx=x1.*mx-abs(x1)-log(1+exp(-2*abs(x1)))+log(2);
elseif (type==4) %%% Laplacian prior - two recified Gaussians
  %Rectified Gaussian
  % Buggy and non-functional...
  t=-x1./sqrt(2*x2);
  erf_table=erfcx(t);
  
  mx=(x1./x2+sqrt(2/pi./x2)./erf_table).*(t<=25)+(-1./x1+2*x2.*x1.^-3-10*x2.^2.*x1.^-5).*(t>25);
  mx2=(x1.^2.*x2.^-2+1./x2+2*x1./x2./sqrt(2*pi*x2)./erf_table).*(t<=25)+(2.*x1.^-2-10*x2.*x1.^-4+74*x2.^2.*x1.^-6).*(t>25);
  Hx=(-log(erfc(min(t,25)))+0.5*log(2*x2/pi)-0.5+x1./sqrt(2*pi*x2)./erf_table).*(t<25)+(log(-x1)-1+2*x2.*x1.^-2-15*x2.^2.*x1.^-4/2+148*x2.^3.*x1.^-6/3).*(t>=25)+0.5*log(pi/2);

  
  
end









