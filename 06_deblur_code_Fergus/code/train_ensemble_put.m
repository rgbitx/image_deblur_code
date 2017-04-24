function x=train_ensemble_put(c,dimensions,x,cx)

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
% Replaces data in the ensemble. Inputs are
%   c          - Class to replace
%   dimensions - The ensemble dimensions matrix
%   x          - Matrix containing all the class data (eg ensemble.mx,
%                ensemble.mx2, etc)
%   cx         - Data for the class
%=========================================================================

if (c>1)
  start=dimensions(1:c-1,1)'*dimensions(1:c-1,2);
else
  start=0;
end
x(start+1:start+dimensions(c,1)*dimensions(c,2))=reshape(cx,1,dimensions(c,1)*dimensions(c,2));
