%--------------------------------------------------------------------------
% File: chooses.m
%
% Goal: Choose the scale value s corresponding to the selected tumor type
%
% Use: [s] = chooses(varargin)
%
% Input: varargin - 1X16 cell containing 8 couples:
%                   'grow_p',value1 - growth law of the primary tumor
%                   'grow_m',value2 - growth law of the matastases
%                   'emission_p',mu_p - colonization coefficient of the
%                                       primary tumor
%                   'emission_m',mu_m - colonization coefficient of the
%                                       metastases
%                   'v_p0', vp0 - volume of the primary tumor at time t=0
%                   'v_m0', vm0 - volume of the newly created metastases
%                   'Vbar', vbar - lower bound of the volume of the 
%                                  metastases whose cumulative number Nv is 
%                                  computed
%                   'tumor_type',value3 - can assume values 'lung' and 'breast'
%
% Output: s - scale value
%
% Comments: These values where chosen in order to have an error less than
%           10^(-6) for a time interval of 400 days. 
%           Notice that for the combination in between two different tumor 
%           growth laws for primary and secondary tumors these values might 
%           not be adequate.
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
function [s] = chooses(varargin)
% default settings
value1 = 'gomp';
value2 = 'gomp';
value3 = 'breast';
control_params = {'grow_p',value1,'grow_m',value2,'tumor_type',value3};
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
switch tumor_type
    case 'lung'
        switch grow_p
            case 'gomp'
                s = 11;
            case 'exp'
                s = 2;
            case 'gen_log'
                s = 4;
            case 'von_ben'
                s = 20;
            case 'power'
                s = 4;
        end
    case 'breast'
        switch grow_p
            case 'gomp'
                s = 11;
            case 'exp'
                s = 5;
            case 'gen_log'
                s = 6;
            case 'von_ben'
                s = 3;
            case 'power'
                s = 4;
        end
end
