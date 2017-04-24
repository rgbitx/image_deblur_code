function c_log_lambda_x=train_ensemble_get_lambda(c,dimensions,log_lambda_x)

% Author: James Miskin
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
% Extracts mixture weight data from the ensemble. Inputs are
%   c            - Class to extract
%   dimensions   - The ensemble dimensions matrix
%   log_lambda_x - Matrix containing all the class data (eg ensemble.log_lambda_x)
%=========================================================================

if (c>1)
  start=sum(prod(dimensions(1:c-1,1:3),2),1);
else
  start=0;
end
c_log_lambda_x=reshape(log_lambda_x(start+1:start+prod(dimensions(c,1:3))),dimensions(c,1),dimensions(c,2),dimensions(c,3));
