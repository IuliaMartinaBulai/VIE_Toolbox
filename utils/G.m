%--------------------------------------------------------------------------
% File: G.m
%
% Goal: Approximate the integral
%             t+gamma
%       G(t) = | h(s,t)ds,  
%              0
%       using the N-point truncated Gauss-Laguerre quadrature rule
%
% Use: [gt] = G(t,h,w,Y,gamma)
%
% Input: t - row array of the evaluation times expressed in days
%        h - integrand of the function G(t)
%        w - weights of the classical Gauss-Laguerre rule
%        Y - matrix ((t(i)+gamma)Exp(-x1(k))_{k=1,...,N, i=1,...,length(t)}
%            with x1 array of N Laguerre zeros
%        gamma - real number. Its value is 0 in the computation of the 
%                total metastatic mass and of the total number of
%                metasteses (vbar=10^-6)
%
% Output: gt - approximation of G(t)
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
function [gt] = G(t,h,w,Y,gamma)
m = length(t);
gt = zeros(1,m);
for i = 1:m
    y = t(i)+gamma;
    r = 1;
    g = w(r)*y*h(Y(r,i),t(i));
    while (r< length(w)) && (abs(g)>0.5e-25)
        gt(i) = gt(i)+g;
        r = r+1;
        g = w(r)*y*h(Y(r,i),t(i));
    end
end
end