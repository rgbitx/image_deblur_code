function [H,f,A,b,L,k,l,u] = qpgen(nx,ne,nc,nb,rank,def,z);

% This function generates a random quadratic program in the form:
%
% min  0.5x'Hx + f'x
%  x
% s.t. Ax  =  b
%      Lx <=  k
%      l  <=  x  <=  u.
%
% Usage is as follows:
%
% [H,f,A,b,L,k,l,u] = qpgen(nx,ne,nc,nb,rank,def,z);
%
% where nx is the dimension of x, ne is the number of equality constraints,
% nc is the number of general inequality constraints and nb is either equal 
% to nx or is zero.
%
% Rank referes to the rank of H and def is 0 for positive semi-definite, 
% and 1 for negative semi-definite. If z is 1 then H will contain zeros.
%
% E.g. copy the following to the Matlab command line to generate a convex
% quadratic programming problem with linear equality, linear inequality
% and simple bound constraints.
%
% n=100; [H,f,A,b,L,k,l,u]=qpgen(n,round(0.5*n),round(1.5*n),n,n,0,0);
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

nx = max(nx,0);
ne = max(ne,0);
nc = max(nc,0);
nb = max(nx,0);

H = []; f = [];
L = []; k = [];
A = []; b = [];
l = []; u = [];

rank = max(rank,0);
rank = min(rank,nx);
if (def > 1 & def < 0), def = 0; end

if z ==1,
  H = [randn(rank,rank);zeros(nx-rank,rank)];
else
  H = randn(nx,rank);
end

if def == 0,
  H = H*H';
else
  H = (H+H')/2;
end
f = 100*randn(nx,1);

if nc,
  L = randn(nc,nx);
  k = rand(nc,1);
end
if ne,
  A = randn(ne,nx);
  b = rand(ne,1);
end
if nb,
  l = -rand(nb,1);
  u = rand(nb,1);
end
