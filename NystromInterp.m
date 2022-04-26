%--------------------------------------------------------------------------
% File: NystromInterp.m
%
% Goal: Compute the Nystrom interpolant, solution of the VIE
%
% Use: [f] = NystromInterp(t,K,g,a,x,w1,E,U)
%
% Input: t - row array of evaluation points
%        K - kernel of the VIE
%        g - right-hand side of the VIE
%        a - solution of the linear system
%        x - Laguerre zeros
%        w1 - weights of the N-point Gauss-Laguerre rule with N=2048
%        E - array exp(-x1) with x1 array of N=2048 Laguerre zeros
%        U - u(x) weight function u computed at the Laguerre zeros x
%
% Output: f - weighted Nystrom interpolant, solution of the VIE
%
% Recalls: ckN.m
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
function [f] = NystromInterp(t,K,g,a,x,w1,E,U)
j = length(a);
Ut = t.^0.25.*exp(-t/2);
Y = t.*E';
% compute the weighted Nystrom interpolant
f = (g(t)+(ckN(t,K,j,x,w1,Y,U)*a)').*Ut;
end