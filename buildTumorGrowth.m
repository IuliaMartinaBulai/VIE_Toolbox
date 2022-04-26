%--------------------------------------------------------------------------
% File: buildTumorGrowth.m
%
% Goal: build the matrix and the right-hand side terms of the linear
%       systems for the metastatic mass and the cumulative number of
%       metastases
%
% Use: [A,b1,b2] = buildTumorGrowth(K,g1,g2,j,x,w1,E,U,gamma)
%
% Input: K - kernel of the VIE            
%        g1 - integrand h(s,t) in the righ-hand side of the VIE 
%             for computing the volume of the metastatic mass
%        g2 - integrand h(s,t) of the righ-hand side of the VIE 
%             for computing the cumulative number of metastases
%        j - size of the solved linear system
%        x - Laguerre zeros
%        w1 - weights of the N-point Gauss-Laguerre rule with N=2048
%        E - array exp(-x1) with x1 array of N=2048 Laguerre zeros
%        gamma - upper limit of integration of the integral at the
%                righ-hand side of the VIE for computing the cumulative
%                number of metastases
%
% Output: A - matrix of the linear system
%         b1 - right-hand side term of the linear system for the volume of 
%              the metastatic mass
%         b2 - right-hand side term of the linear system for the cumulative
%              number of metastases
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
function [A,b1,b2] = buildTumorGrowth(K,g1,g2,j,x,w1,E,U,gamma)
Y = x(1:j).*E';
Y1 = (x(1:j)+gamma).*E';
% build the matrix A of the linear system
A = eye(j)-ckN(x(1:j),K,j,x,w1,Y,U).*U';
% build the right-hand side term b of the linear system for the volume of 
% the metastatic mass
b1 = (G(x(1:j),g1,w1,Y,0).*U)';
% build the right-hand side term b of the linear system for the cumulative 
% number of metastases
b2 = (G(x(1:j),g2,w1,Y1,gamma).*U)';
end