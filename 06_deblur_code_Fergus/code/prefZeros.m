function [ res ] = prefZeros( n, z )

%
%  [ res ] = prefZeros( n, z )
%
%  res is a string of length z, containing the number n plus leading
%  zeros.
%

%
%  (c) in 1999 by Markus Weber
%
%  Califonia Institute of Technology
% 

%
%  $Id: prefZeros.m,v 1.1.1.1 2003/10/20 16:22:42 fergus Exp $
%
%  $Log: prefZeros.m,v $
%  Revision 1.1.1.1  2003/10/20 16:22:42  fergus
%  ECCV '04 code
%
%  Revision 1.1.1.1  2002/11/21 19:01:47  fergus
%  Inital import of constellation model code
%



res = fliplr(['000000000000' int2str(n)]);

res = fliplr(res(1:z));
