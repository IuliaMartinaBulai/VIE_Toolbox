%--------------------------------------------------------------------------
% File: ckN.m
%
% Goal: Approximate the integrals 
%               +oo
%    c_k(t) = t | K(t exp(-z),t)l_{n+1,k}(t exp(-z)) exp(-z)dz,  k=1,...,j
%               0 
%       using the N-point truncated Gauss-Laguerre quadrature rule
%
% Use: [c] = ckN(t,K,j,x,w1,Y,U)
% 
% Input: t - row array of evaluation points
%        K - kernel of the VIE
%        j - size of the solved linear system
%        x - Laguerre zeros
%        w1 - weights of the N-point Gauss-Laguerre rule 
%        Y - matrix t.*exp(-x1)' with x1 array of N Laguerre zeros
%        U - u(x) weight function u computed at the Laguerre zeros x
%
% Output: c - matrix (c_k(t_i))_{i=1,...,length(t), k=1,...,j}
%
% Recalls: lnk.m
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
function [c] = ckN(t,K,j,x,w1,Y,U)
n = length(x);
x = [x,4*n];
m = length(t);
c = zeros(m,j);
for i = 1:m
    for k = 1:j
        % compute the approximation of c_k(t_i)
        nu = 1;
        cn = t(i)*w1(nu)*K(Y(nu,i),t(i))*lnk(k,x,Y(nu,i))/U(k);
        while (nu < length(w1)) && (abs(cn) > 0.5e-25)
            c(i,k) = c(i,k)+cn;
            nu = nu+1;
            cn = t(i)*w1(nu)*K(Y(nu,i),t(i))*lnk(k,x,Y(nu,i))/U(k);
        end
    end
end
end