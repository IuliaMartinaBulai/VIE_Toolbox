%--------------------------------------------------------------------------
% File: NystromMethodTumorGrowth.m
%
% Goal: Compute the volume of the metastatic mass and the cumulative
%       number of metastases using the Nystrom method in [BDBLS]
%
% Use: [M,Nv,j,C] = NystromMethodTumorGrowth(t,g,n,gamma,theta)
%
% Input: t - row array of the evaluation times expressed in days
%        g - 1X4 cell with
%            g{1} integrand h(s,t) in the right-hand side of the VIE 
%            for computing the volume of the metastatic mass
%            g{2} integrand h(s,t) of the right-hand side of the VIE 
%            for computing the cumulative number of metastases
%            g{3} kernel K(s,t) of the VIEs
%            g{4} weight function u(s) of the weighted space 
%        n - number of knots
%        gamma - upper limit of integration of the integral at the
%                right-hand side of the VIE for computing the cumulative
%                number of metastases
%        theta - truncation parameter
%
% Output: M - volume of the total metastatic burden at t days in mm^3
%         Nv - cumulative number of metastases whose volume is larger than 
%              Vbar at t days
%         j - size of the solved linear system
%         C - condition number of the solved linear system
%
% Recalls: gaussq.m, buildTumorGrowth.m, NystromInterpTumorGrowth.m
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
function [M,Nv,j,C] = NystromMethodTumorGrowth(t,g,n,gamma,theta)
% compute n Laguerre knots 
[x,~] = gaussq(6,n,0,0,0,[0,0]);
% compute N=2048 Laguerre knots and weights
[x1,w1] = gaussq(6,2048,0,0,0,[0,0]);
E1 = exp(-x1);
j = 1;
while (j <= n) && (x(j) < 4*n*theta)
    j = j+1;
end
j = j-1;
g1 = g{1};
g2 = g{2};
K = g{3};
u = g{4};
U = u(x(1:j));
% build the matrix and the right-hand side terms of the linear systems for
% the total metastatic mass and the cumulative number of metastases
[A,b1,b2] = buildTumorGrowth(K,g1,g2,j,x,w1,E1,U,gamma);
% compute the condition number of the matrix A
C = cond(A,inf);
% compute the solution of the linear systems 
a1 = A\b1;
a2 = A\b2;
% compute the volume of the metastatic mass and the cumulative number of
% metastases
[M,Nv] = NystromInterpTumorGrowth(t,K,g1,g2,a1,a2,x,w1,E1,U,gamma);
end