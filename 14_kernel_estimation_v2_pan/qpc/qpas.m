%  This dual active-set algorithm attempts to solve
%  quadratic programming problems in the following form:
%
%   x* = arg min   0.5x'Hx + f'x
%             x
%
%            s.t.  Ax  = b,        linear equality constraints,
%                  Lx <= k,        general linear inequality constraints,
%                  l  <= x <= u,   bound constraints.
%
%  In Matlab a call to this function would look like:
%
%    [x,err,lm] = qp(H,f,L,k,A,b,l,u,display);
%
%  where the inputs are:
%
%        H: An (n x n) positive semi-definite symmetric matrix
%
%	     f: A n element column vector
%	
%    (L,k): General linear inequality constraints
%
%    (A,b): General linear equality constraints
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
%    lm.equality: Lagrange multipliers for equality constraints
%
%  lm.inequality: Lagrange multipliers for general inequality constraints
%
%  lm.lowerbound: Lagrange multipliers for lower bound constraints
%
%  lm.upperbound: Lagrange multipliers for upper bound constraints
%
%  If (A,b) and/or (L,k) and/or l and/or u, are not being used, then
%  set the respective entry to [] (i.e. the empty matrix).
%
%  E.g.1 suppose that only L and k are being used then a call would look
%  like
%
%    x = qp(H,f,L,k);
%
%  E.g.2 suppose that (A,b) and l are required, then a call would look
%  like
%
%    x = qp(H,f,[],[],A,b,l);
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