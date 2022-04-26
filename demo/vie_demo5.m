%--------------------------------------------------------------------------
% vie_demo5:  
%
% This demo solves the VIE in order to compute the metastatic mass, M(t), 
% and the cumulative number of metastases, Nv(t), for breast tumor data.
% We consider and compare five different tumor growth laws to describe the 
% tumor metastatic growth: exponential, power-law, Gompertz, generalized 
% logistic and von Bertalanffy-West laws.
%
% This reproduces figure 2 from "Modeling metastatic tumor evolution, 
% numerical resolution and growth prediction", IM Bulai, MC De Bonis, C
% Laurita, V Sagaria, 2022.
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
%% Load numerical results for breast, five different growth laws
clc
clear all
close all

fprintf('Welcome to VIE demo #5\n');
addpath(fullfile('','Data'));
fprintf('Loading the numerical results for breast \n');
load('mass_number_breast_t_60.mat');

t = [1:T];

fprintf('Plot metastatic mass \n');
figure
semilogy(t,M_1, 'r', t,M_2,'b', t,M_3,'g', t,M_4,'c', t,M_5,'k','linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
ylabel('Metastatic mass (volume)')
legend('gomp', 'exp', 'gen log', 'von bert', 'power law', 'Location','best')

pause

figure
fprintf('Plot cumulative number of metastases \n');
semilogy(t,N_1, 'r', t,N_2,'b', t,N_3,'g', t,N_4,'c', t,N_5,'k', 'linewidth',2)
set(gca,'fontsize',16)
title('Breast')
xlabel('Time (days)')
ylabel('Cumulative number')
legend('gomp', 'exp', 'gen log', 'von bert', 'power law', 'Location','best')