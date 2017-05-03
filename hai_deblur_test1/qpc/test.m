% Run this file to test routines
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


% Generate a strictly convex quadratic program (random)
n=40; [H,f,A,b,L,k,l,u]=qpgen(n,round(0.2*n),3*n,n,n,0,0);


% clear all
% 
% load test_data

dsp=0;

% Test qp
disp(sprintf('\n'));
disp('Testing qpas routine');
try, 
    tic;xqpas=qpas(H,f,L,k,A,b,l,u,dsp);toc
    disp('Test OK');
catch, 
    disp('There was an error'); 
end

% Test qpip
disp(sprintf('\n'));
disp('Testing qpip routine');
try, 
    tic;xqpip=qpip(H,f,L,k,A,b,l,u,dsp); toc
    disp('Test OK');
catch, 
    disp('There was an error'); 
end

% Test qps_as
disp(sprintf('\n'));
disp('Testing qps_as routine');
try, 
    tic;xqps_as=qps_as(H,f,l,u,dsp);toc
    disp('Test OK');
catch, 
    disp('There was an error'); 
end

% Test qps_ip
disp(sprintf('\n'));
disp('Testing qps_ip routine');
try, 
    tic;xqps_ip=qps_ip(H,f,l,u,dsp);toc
    disp('Test OK');
catch, 
    disp('There was an error'); 
end

% Test qps_mq
disp(sprintf('\n'));
disp('Testing qps_mq routine');
try, 
    tic;xqps_mq=qps_mq(H,f,l,u);toc
    disp('Test OK');
catch, 
    disp('There was an error'); 
end


disp(sprintf('\n'));