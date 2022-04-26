%--------------------------------------------------------------------------
% File: NystromInterpTumorGrowth.m
%
% Goal: Compute the Nystrom interpolants, solutions of VIEs, approximating
%       the volume of the metastatic mass and the cumulative number of
%       metastases
%
% Use: [M,Nv] =  NystromInterpTumorGrowth(t,K,g1,g2,a1,a2,x,w1,E,U,gamma)
%
% Input: t -  row array of the evaluation times expressed in days
%        K - kernel of the VIE
%        g1 - integrand h(s,t) in the righ-hand side of the VIE 
%             for computing the volume of the metastatic mass
%        g2 - integrand h(s,t) of the righ-hand side of the VIE 
%             for computing the cumulative number of metastases
%        a1 - solution of the linear system for computing the volume of the
%             metastatic mass
%        a2 - solution of the linear system for computing the cumulative 
%             number of metastases
%        x - Laguerre zeros
%        w1 - weights of the N-point Gauss-Laguerre rule with N=2048
%        E - array exp(-x1) with x1 array of N=2048 Laguerre zeros
%        U - u(x) weight function u computed ad the Laguerre zeros x
%        gamma - upper limit of integration of the integral at the
%                righ-hand side of the VIE for computing the cumulative
%                number of metastases
%
% Output: M - Nystrom interpolant, solution of the VIE for computing
%             the volume of the metastatic mass
%         Nv - Nystrom interpolant, solution of the VIE for computing
%              the cumulative number of metastases
%
% Recalls: ckN.m, G.m
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
function [M,Nv] = NystromInterpTumorGrowth(t,K,g1,g2,a1,a2,x,w1,E,U,gamma)
j=length(a1);
Y = t.*E';
Y1 = (t+gamma).*E';
% compute the modified moments using a N-point truncated Gauss-Laguerre 
% quadrature rule
C = ckN(t,K,j,x,w1,Y,U);
% compute the Nystrom interpolant for the volume of the metastatic mass 
M = G(t,g1,w1,Y,0)'+C*a1; 
% compute the Nystrom interpolant for the cumulative number of metastases
Nv = G(t,g2,w1,Y1,gamma)'+C*a2;
end