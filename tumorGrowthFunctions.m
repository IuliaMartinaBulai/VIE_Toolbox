%--------------------------------------------------------------------------
% File: tumorGrowthFunctions.m
%
% Goal: Create a database of possible tumor growth laws both for primary
%       and metastatic growth, for lung and breast case studies
%
% Use: [s,gamma,g] = tumorGrowthFunctions(varargin)
%
% Input: varargin - 1X16 cell containing 8 couples:
%                   'grow_p',value1 - growth law of the primary tumor
%                   'grow_m',value2 - growth law of the metastases
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
% Output: g - 1X4 cell with
%             g{1} iintegrand h(s,t) in the right-hand side of the VIE 
%             for computing the volume of the metastatic mass
%             g{2} integrand h(s,t) of the right-hand side of the VIE 
%             for computing the cumulative number of metastases
%             g{3} kernel K(s,t) of the VIEs
%             g{4} weight function u(s) of the weighted space 
%         gamma - upper limit of integration of the integral at the
%                right-hand side of the VIE for computing the cumulative
%                number of metastases
%         s - scale parameter of the time t
%
% Recalls: chooses.m
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
function [g,gamma,s] = tumorGrowthFunctions(varargin)
% default growth law for the primary tumor
value1 = 'gomp';
% default growth law for the metastases
value2 = 'gomp';
% default case study
value3 = 'breast';
% default parameter values
al = 2/3;
mu_p = 10^-3; % day^(-1)mm^(-3 al)
mu_m = 10^-3; % day^(-1)mm^(-3 al)
vp0 = 1;  % mm^3  
vm0 = 10^-6; % mm^3  
vbar = 10^-6; % mm^3
control_params = {'grow_p',value1,'grow_m',value2, ...
    'emission_p',mu_p,'emission_m',mu_m,...
    'v_p0', vp0,'v_m0', vm0,'Vbar', vbar,'tumor_type',value3};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

[s] = chooses('grow_p',value1,'grow_m',value2,'tumor_type',value3);

switch tumor_type
    case 'lung'
        switch grow_p
            case 'gomp'
                a = 0.743*s;
                bet = 0.0792*s;
                vp = @(t) (v_p0*exp(a*(1-exp(-bet*t))/bet));
            case 'exp'
                vp0 = 13.2;
                a0 = 0.257*s;
                vp = @(t) (vp0*exp(a0*t));
            case 'gen_log'
                a = 2555*s;
                K = 4378;
                nu = 1.4e-04;
                vp = @(t) (v_p0*K/(v_p0^nu+(K^nu-v_p0^nu)*exp(-a*nu*t))^(1/nu));
            case 'von_ben'
                a = 7.72*s;
                gamm = 0.947;
                b = 6.75*s;
                vp = @(t) ((a/b+(v_p0^(1-gamm)-a/b)*exp(-b*(1-gamm)*t))^(1/(1-gamm)));
            case 'power'
                a = 0.921*s;
                gamm = 0.788;
                vp = @(t) ((v_p0^(1-gamm)+(1-gamm)*a*t)^(1/(1-gamm)));
        end
        switch grow_m
            case 'gomp'
                a = 0.743*s;
                bet = 0.0792*s;
                vm = @(t) (v_m0*exp(a*(1-exp(-bet*t))/bet));
                if Vbar > v_m0*exp(a/bet) || Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = log (1-bet/a*log(Vbar/v_m0))/bet;
            case 'exp'
                a0 = 0.257*s;
                vm = @(t) (v_m0*exp(a0*t));
                if Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = -1/a0*log(Vbar/v_m0);
            case 'gen_log'
                a = 2555*s;
                K = 4378;
                nu = 1.4e-04;
                vm = @(t) (v_m0*K/(v_m0^nu+(K^nu-v_m0^nu)*exp(-a*nu*t))^(1/nu));
                if (Vbar > K || v_m0>K)  || Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = log ((v_m0^nu*((K/Vbar)^nu-1))/(K^nu-v_m0^nu))/(a*nu);
            case 'von_ben'
                a = 7.72*s;
                gamm = 0.947;
                b = 6.75*s;
                vm = @(t) ((a/b+(v_m0^(1-gamm)-a/b)*exp(-b*(1-gamm)*t))^(1/(1-gamm)));
                if (Vbar > (a/b)^(1/(1-gamm)) || v_m0 > (a/b)^(1/(1-gamm)) )  || Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = log ((a/b-Vbar^(1-gamm))/(a/b-v_m0^(1-gamm)))/(b*(1-gamm));
            case 'power'
                a = 0.921*s;
                gamm = 0.788;
                vm = @(t) ((v_m0^(1-gamm)+(1-gamm)*a*t)^(1/(1-gamm)));
                if Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = (v_m0^(1-gamm)-Vbar^(1-gamm))/((1-gamm)*a);
        end
    case 'breast'
        switch grow_p
            case 'gomp'
                a = 0.56*s;
                bet = 0.0719*s;
                vp = @(t) (v_p0*exp(a*(1-exp(-bet*t))/bet));
            case 'exp'
                vp0 = 68.2;
                a0 = 0.0846*s;
                vp = @(t) (vp0*exp(a0*t));
            case 'gen_log'
                a = 2753*s;
                K = 1964;
                nu = 2.68e-05;
                vp = @(t) (v_p0*K/(v_p0^nu+(K^nu-v_p0^nu)*exp(-a*nu*t))^(1/nu));
            case 'von_ben'
                a = 2.32*s;
                gamm = 0.918;
                b = 0.808*s;
                vp = @(t) ((a/b+(v_p0^(1-gamm)-a/b)*exp(-b*(1-gamm)*t))^(1/(1-gamm)));
            case 'power'
                a = 1.32*s;
                gamm = 0.58;
                vp = @(t) ((v_p0^(1-gamm)+(1-gamm)*a*t)^(1/(1-gamm)));
        end
        switch grow_m
            case 'gomp'
                a = 0.56*s;
                bet = 0.0719*s;
                vm = @(t) (v_m0*exp(a*(1-exp(-bet*t))/bet));
                if Vbar > v_m0*exp(a/bet) || Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = log (1-bet/a*log(Vbar/v_m0))/bet;
            case 'exp'
                a0 = 0.0846*s;
                vm = @(t) (v_m0*exp(a0*t));
                if Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = -1/a0*log(Vbar/v_m0);
            case 'gen_log'
                a = 2753*s;
                K = 1964;
                nu = 2.68e-05;
                vm = @(t) (v_m0*K/(v_m0^nu+(K^nu-v_m0^nu)*exp(-a*nu*t))^(1/nu));
                if (Vbar > K || v_m0>K)  || Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = log ((v_m0^nu*((K/Vbar)^nu-1))/(K^nu-v_m0^nu))/(a*nu);
            case 'von_ben'
                a = 2.32*s;
                gamm = 0.918;
                b = 0.808*s;
                vm = @(t) ((a/b+(v_m0^(1-gamm)-a/b)*exp(-b*(1-gamm)*t))^(1/(1-gamm)));
                if (Vbar > (a/b)^(1/(1-gamm)) || v_m0 > (a/b)^(1/(1-gamm)) )  || Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = log ((a/b-Vbar^(1-gamm))/(a/b-v_m0^(1-gamm)))/(b*(1-gamm));
            case 'power'
                a = 1.32*s;
                gamm = 0.58;
                vm = @(t) ((v_m0^(1-gamm)+(1-gamm)*a*t)^(1/(1-gamm)));
                if Vbar < v_m0
                    error('V_{bar} not feasible ');
                end
                gamma = (v_m0^(1-gamm)-Vbar^(1-gamm))/((1-gamm)*a);
        end
end
% emission rate function for the primary tumor
be_p = @(t) s*emission_p*t^al;
% emission rate function for the metastases
be_m = @(t) s*emission_m*t^al;
% h(s,t) in the paper, for the volume of the metastatic mass
g{1} = @(s,t) vm(t-s)*be_p(vp(s)); 
% h(s,t) in the paper, for the cumulative number of metasteses
g{2} = @(s,t) be_p(vp(s));         
% k(s,t) in the paper, the kernel of the VIE
g{3} = @(s,t) be_m(vm(t-s));   
% u(s) in the paper, the weight function of the weighted space where
% the solutions of the VIEs lives
g{4} = @(y) y.^(1/4).*exp(-y/2);   
end