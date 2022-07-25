%--------------------------------------------------------------------------
% vie_demo1:  
%
% This demo solves the VIE in order to compute the metastatic mass, M(t), 
% and the cumulative number of metastases, Nv(t), for breast tumor data
% assuming that both the primary and secondary tumors growth laws are 
% Gompertz laws. 
%
% This reproduces figure 2 from "MatLab Toolbox for the numerical solution
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
fprintf('Welcome to VIE demo #1\n');
addpath('..');
addpath('../Utils');

fprintf('Solve the VIE for default values \n');
fprintf('Case study: breast \n');
T = 60;
tic
[M,err_rel_M,Nv,err_rel_Nv,j,C] = tumorGrowth(T);

fprintf(' Index J and condition number C \n');
[j, C]
errM = max(err_rel_M);
errNv = max(err_rel_Nv);
fprintf('Error for metastatic mass and for cumulative number of metastases \n');

t = [1:T];

fprintf('Plot the error for the metastatic mass \n');
figure
plot(t,err_rel_M, 'linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
legend('Error metastatic mass','Location','Best')

pause 

fprintf('Plot the error for the cumulative number of metastases \n');
figure
plot(t,err_rel_Nv, 'linewidth',2);
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
legend('Error cumulative number','Location','Best')

pause 

fprintf('Plot metastatic mass \n');
figure
plot(t,M, 'linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')

pause

fprintf('Plot cumulative number of metastases \n');
figure
plot(t,Nv, 'linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
ylabel('Cumulative number')
toc
