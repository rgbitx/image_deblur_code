%  This primal-dual predictor-corrector algorithm attempts to solve
%  quadratic programming problems in the following form:
%
%   x* = arg min   0.5x'Hx + f'x
%             x
%
%            s.t.  l  <= x <= u,   bound constraints.
%
%  In Matlab a call to this function would look like:
%
%    [x,err,lm] = qps_ip(H,f,l,u,display);
%
%  where the inputs are:
%
%        H: An (n x n) positive semi-definite symmetric matrix
%
%	     f: A n element column vector
%	
%        l: Element-wise lower bound constraints
%
%        u: Element-wise upper bound constraints
%
%  display: If display>0 then iteration information is displayed
%
%
%  and the outputs are:
%
%              x: the optimal solution (if obtained)
%
%            err: error number, if err=0, then x is optimal
%
%             lm: structure of Lagrange multipliers 
% 
%  lm.lowerbound: Lagrange multipliers for lower bound constraints
%
%  lm.upperbound: Lagrange multipliers for upper bound constraints
%
% NOTE: Both l and u must be provided, however, +/-inf entries are
%       acceptable so that "free" variables are catered for.
%  
%   +----------------------------------------------+
%   | Written by Adrian Wills,                     |
%   |            School of Elec. Eng. & Comp. Sci. |
%   |            University of Newcastle,          |
%   |            Callaghan, NSW, 2308, AUSTRALIA   |
%   |                                              |
%   | Last Revised  25 May 2007.                   |
%   |                                              |
%   | Copyright (C) Adrian Wills.                  |
%   +----------------------------------------------+
%  
% The current version of this software is free of charge and 
% openly distributed, BUT PLEASE NOTE:
% 
% This software must be referenced when used in a published work.
% 
% This software may not be re-distributed as a part of a commercial product. 
% If you distribute it in a non-commercial products, please contact me first, 
% to make sure you ship the most recent version.
% 
% This software is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% 
% IF IT FAILS TO WORK, IT'S YOUR LOSS AND YOUR PROBLEM.