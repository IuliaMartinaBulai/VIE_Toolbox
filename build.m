%--------------------------------------------------------------------------
% File: build.m
%
% Goal: Build the matrix and the right-hand side term of the linear system 
%
% Use: [A,b] = build(K,g,j,x,w1,E,U)
%
% Input: K - kernel of the VIE
%        g - right-hand side of the VIE
%        j - size of the solved linear system
%        x - Laguerre zeros
%        w1 - weights of the N-point Gauss-Laguerre rule with N=2048
%        E - array exp(-x1) with x1 array of N=2048 Laguerre zeros
%        U - u(x) weight function u computed at the Laguerre zeros x
%
% Output: A - matrix of the linear system
%         b - right-hand side term of the linear system
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
function [A,b] = build(K,g,j,x,w1,E,U)
Y = x(1:j).*E';
% build the matrix A of the linear system
A = eye(j)-ckN(x(1:j),K,j,x,w1,Y,U).*U';
% build the right-hand side term b of the linear system
b = (g(x(1:j)).*U)';
end