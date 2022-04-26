%--------------------------------------------------------------------------
% File: vie.m
%
% Goal: Compute the solution of the following Volterra Integral Equation (VIE) 
% of the second kind on infinite intervals 
%             t
%      f(t) - | K(s,t)f(s)ds=g(t),    t>0,
%             0 
% and the corresponding weighted approximation error
%
% Use: [f,err,j,C] = vie(t,K,g)
%
% Input: t - row array of evaluation points
%        K - kernel of the VIE
%        g - right-hand side of the VIE
%
% Output: f - solution of the VIE at t
%         err - weighted absolute error
%         j - size of the solved linear system
%         C - condition number of the solved linear system
%
% Recalls: NystromMethodVIE.m
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
function [f,err,j,C] = vie(t,K,g)
u = t.^(1/4).*exp(-t/2);
theta = 0.25;
k = 9;
[f_0,~,~] = NystromMethodVIE(t,K,g,2^k,theta);
[f,j,C] = NystromMethodVIE(t,K,g,2^(k-1),theta);
% control of the condition number of the solved linear system 
if C > 1/eps
    sprintf('Warning: Results may be inaccurate. COND = %d',C)
end
% compute the weigthed absolute error
err = max(abs(f_0*u'-f*u'));
% control of the convergence of the method using the weigthed absolute error
if err > 1
    disp('The method is not convergent');
    f = [];
end
end