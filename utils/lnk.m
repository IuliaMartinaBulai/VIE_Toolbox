%--------------------------------------------------------------------------
% File: lnk.m
%
% Goal: Compute the k-th fundamental Lagrange polynomial based on the 
%       interpolation knots x at the point t
%
% Use: [l] = lnk(k,x,t)
%
% Input: k - index of the fundamental Lagrange polynomial
%        x - interpolation knots
%        t - evaluation point
%
% Output: l - k-th fundamental Lagrange polynomial at t
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
function [l] = lnk(k,x,t)
n = length(x);
% compute the k-th fundamental Lagrange polynomial 
l = prod([(t-x(1:k-1))./(x(k)-x(1:k-1)),(t-x(k+1:n))./(x(k)-x(k+1:n))]);
end