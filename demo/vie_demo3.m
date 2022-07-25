%--------------------------------------------------------------------------
% vie_demo3:  
%
% This demo computes the weighted solution of a VIE
%
% This reproduces figure 4 from "MatLab Toolbox for the numerical solution
% of linear Volterra integral equations arising in metastatic tumor growth
% models", IM Bulai, MC De Bonis, C Laurita, V Sagaria, 2022. 
%
% Authors: IM Bulai, MC De Bonis, C Laurita, V Sagaria
% Date last modified: April, 2022
%
% This file is part of the VIE toolbox (Volterra Integral Equation toolbox)
% Copyright (C) 2022, IM Bulai, MC De Bonis, C Laurita, V Sagaria.
%
% The VIE toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation.
%
% The VIE toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the VIE toolbox.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
%% 
clc
clear all
close all

fprintf('Welcome to VIE demo #3\n');
addpath('..');
addpath('../Utils');
fprintf('Solve the VIE  \n');

% kernel, right-hand side term and Laguerre weight
K=@(s,t) s./(s.^2 + t.^2 + 1);
g=@(t) sqrt(1+t.^2).*atan(t./sqrt(1+t.^2));
u=@(y) y.^(1/4).*exp(-y/2);
% exact solution 
f=@(t) t;

% row array of the evaluation points
t=linspace(0,40);

tic
[fm,err,j,C]=vie(t,K,g);

fprintf('Index J, condition number and error\n');
[j,C,err]

fprintf('Plot the Nystrom interpolating function \n');
figure
plot(t,fm, 'linewidth',2)
set(gca,'fontsize',16)
title('Weighted Nystrom interpolating function')
xlabel('t')
ylabel('u(t)F_{n,N}(t)')
toc
