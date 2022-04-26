%--------------------------------------------------------------------------
% File: NystromMethodVIE.m
%
% Goal: Compute the solution of the VIE using the Nystrom method in [BDBLS]
%
% Use: [f,j,C] = NystromMethodVIE(t,K,g,n,theta)
%
% Input:  t - row array of evaluation points
%         K - kernel of the VIE 
%         g - right-hand side of the VIE         
%         n - number of knots
%         theta - truncation parameter
%       
%
% Output: f - solution of the VIE 
%         j - size of the solved linear system
%         C - condition number of the solved linear system
%
% Recalls: gaussq.m, build.m, NystromInterp.m
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
function [f,j,C] = NystromMethodVIE(t,K,g,n,theta)
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
U = x(1:j).^0.25.*exp(-x(1:j)/2);
% build the matrix and the right-hand side term of the linear system
[A,b] = build(K,g,j,x,w1,E1,U);
% compute the condition number of the matrix A
C = cond(A,inf);
% compute the solution of the linear system
a = A\b;
% compute the solution of the VIE
f = NystromInterp(t,K,g,a,x,w1,E1,U);
end