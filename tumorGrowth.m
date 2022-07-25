%--------------------------------------------------------------------------
% File: tumorGrowth.m
%
% Goal: Compute the volume of the metastatic mass and the cumulative number 
%       of metastases as solutions of Volterra Integral Equations
%
% Use: [M,err_rel_M,Nv,err_rel_Nv,j,C] = tumorGrowth(T,varargin{:})
%
% Input: T - final time point expressed in days
%        varargin - 1X16 cell containing 8 couples:
%                   'grow_p',value1 - growth law of the primary tumor
%                   'grow_m',value2 - growth law of the metastases
%                   'emission_p',mu_p - colonization coefficient of the
%                                       primary tumor
%                   'emission_m',mu_m  - colonization coefficient of the
%                                        metastases
%                   'v_p0', vp0 - volume of the primary tumor at time t=0
%                   'v_m0', vm0 - volume of the newly created metastases
%                   'Vbar', vbar - lower bound of the volume of the 
%                                  metastases whose cumulative number Nv is 
%                                  computed
%                   'tumor_type', value3 - can assume values 'lung' and 'breast'
%
% Output: M - 1xT array of the approximate values of the metastatic mass 
%             at t=[1:T] days
%         err_rel_M - 1xT array of the relative errors related to the 
%                     approximate values in M
%         Nv - 1xT array of the approximate values of the cumulative number 
%              of metastases whose volume is larger than Vbar at t=[1:T] days
%         err_rel_Nv - 1xT array of the relative errors related to the
%                     approximate values in Nv
%         j - dimension of the solved linear system
%         C - condition number of the solved linear system
%
% Recalls: tumorGrowthFunctions.m, NystromMethodTumorGrowth.m
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
function [M,err_rel_M,Nv,err_rel_Nv,j,C] = tumorGrowth(T,varargin)
[g,gamma,s] = tumorGrowthFunctions(varargin{:});
t = (1:T)/s;
k = 9;
theta = 0.25;
[M_0,Nv_0,~,~] = NystromMethodTumorGrowth(t,g,2^k,gamma,theta);
[M,Nv,j,C] = NystromMethodTumorGrowth(t,g,2^(k-1),gamma,theta);
% control of the condition number of the solved linear system 
if C > 1/eps
    sprintf('Warning: Results may be inaccurate. COND = %d',C)
end
% compute the relative errors related to the approximate values of the 
% volume of the metastatic mass
err_rel_M = abs(M_0-M)./abs(M_0);
errM = max(err_rel_M);
% compute the relative errors related to the approximate values of the 
% cumulative number of metasteses 
err_rel_Nv = abs(Nv_0-Nv)./abs(Nv_0);
errNv = max(err_rel_Nv);
% control of the convergence of the method using the weigthed absolute
% errors
flagM = 0;
if errM > 1
    flagM = 1;
end
flagNv = 0;
if errNv > 1
    flagNv = 1;
end
if flagM && flagNv
    disp('The method is not convergent');
    M = [];
    Nv = [];
elseif flagM
    disp('The method is not convergent for the approximation of M');
    M = [];
elseif flagNv
    disp('The method is not convergent for the approximation of N');
    Nv = [];
end
end
