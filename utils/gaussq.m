%--------------------------------------------------------------------------
% File: gaussq.m
%
%        This set of routines computes the nodes t(j) and weights
%        w(j) for gaussian-type quadrature rules with pre-assigned
%        nodes.  These are used when one wishes to approximate
%
%                 integral (from a to b)  f(x) w(x) dx
%
%                              n
%        by                   sum w  f(t )
%                             j=1  j    j
%
%        (note w(x) and w(j) have no connection with each other.)
%        here w(x) is one of six possible non-negative weight
%        functions (listed below), and f(x) is the
%        function to be integrated.  Gaussian quadrature is particularly
%        useful on infinite intervals (with appropriate weight
%        functions), since then other techniques often fail.
%
%        Associated with each weight function w(x) is a set of
%        orthogonal polynomials.  The nodes t(j) are just the zeros
%        of the proper n-th degree polynomial.
%
%     input parameters (all real numbers are in double precision)
%
%        kind     an integer between 1 and 6 giving the type of
%                 quadrature rule:
%
%        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
%        kind = 2:  chebyshev quadrature of the first kind
%                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
%        kind = 3:  chebyshev quadrature of the second kind
%                   w(x) = sqrt(1 - x*x) on (-1, 1)
%        kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
%                   (-infinity, +infinity)
%        kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
%                   beta on (-1, 1), alpha, beta .gt. -1.
%                   note: kind=2 and 3 are a special case of this.
%        kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*
%                   x**alpha on (0, +infinity), alpha .gt. -1
%
%        n        the number of points used for the quadrature rule
%        alpha    real parameter used only for gauss-jacobi and gauss-
%                 laguerre quadrature (otherwise use 0.d0).
%        beta     real parameter used only for gauss-jacobi quadrature--
%                 (otherwise use 0.d0)
%        kpts     (integer) normally 0, unless the left or right end-
%                 point (or both) of the interval is required to be a
%                 node (this is called gauss-radau or gauss-lobatto
%                 quadrature).  then kpts is the number of fixed
%                 endpoints (1 or 2).
%        endpts   real array of length 2.  contains the values of
%                 any fixed endpoints, if kpts = 1 or 2.
%
%     output parameters (both double precision arrays of length n)
%
%        t        will contain the desired nodes.
%        w        will contain the desired weights w(j).
%
%     underflow may sometimes occur, but is harmless.
%
%     references
%        1.  golub, g. h., and welsch, j. h., "calculation of gaussian
%            quadrature rules," mathematics of computation 23 (april,
%            1969), pp. 221-230.
%        2.  golub, g. h., "some modified matrix eigenvalue problems,"
%            siam review 15 (april, 1973), pp. 318-334 (section 7).
%        3.  stroud and secrest, gaussian quadrature formulas, prentice-
%            hall, englewood cliffs, n.j., 1966.
%
%        original version 20 jan 1975 from stanford
%        modified 21 dec 1983 by eric grosse
%          imtql2 => gausq2
%          hex constant => d1mach (from core library)
%          compute pi using datan
%          removed accuracy claims, description of method
%          added single precision version
%
%        Translation of the fortran procedure
%        http://www.netlib.org/go/gaussq.f
%--------------------------------------------------------------------------
function [t,w]=gaussq(kind,n,alpha,beta,kpts,endpts)
[a1,b1,mu1]=class1(kind,n,alpha,beta);
if kpts == 1
    [gam]=solve(endpts(1),n,a1,b1);
    a1(n)=gam*b1(n-1)^2+endpts(1);
elseif kpts==2
    [gam]=solve(endpts(1),n,a1,b1);
    [gam1]=solve(endpts(2),n,a1,b1);
    t1=((endpts(1)-endpts(2))/(gam1-gam));
    b1(n-1)=sqrt(t1);
    a1(n)=gam*t1+endpts(1);
end
[t,~,w,~]=gausq2(n,a1,b1);
w=w.^2*mu1;

function[d,e,z,ierr]= gausq2(n, d, e)
% This subroutine is a translation of an algol procedure, num. math.12,
% 377-383(1968) by Martin and Wilkinson, as modified in num. math. 15,
% 450(1970) by Dubrulle. Handbook for auto. comp., vol.ii-linear algebra,
% 241-248(1971). This is a modified version of the 'eispack' routine imtql2.
% This subroutine finds the eigenvalues and first components of the
% eigenvectors of a symmetric tridiagonal matrix by the implicit ql method.
%
% Input: n - the order of the matrix
%        d - the diagonal elements of the input matrix
%        e - the subdiagonal elements of the input matrix
%            in its first n-1 positions. e(n) is arbitrary;
%        z - the first row of the identity matrix
%
% Output: d - the eigenvalues in ascending order, if an error exit is made,
%             the eigenvalues are correct but unordered for indices
%             1, 2, ..., ierr-1;
%         e - has been destroyed;
%         z - the first components of the orthonormal eigenvectors of the
%             symmetric tridiagonal matrix. If an error exit is made,
%             z contains the eigenvectors associated with the stored eigenvalues;
%         ierr - is set to zero for normal return, j if the j-th eigenvalue
%                has not been determined after 30 iterations.
%     ------------------------------------------------------------------
z=zeros(1,n); z(1)=1;
ierr = 0;
if  n ~= 1
    e(n) = 0;
    for l = 1 : n
        j = 0;
        % look for small sub-diagonal element
        ivar =1;
        while ivar ==1
            for  m = l : n
                if m == n
                    break
                end
                if abs(e(m)) <= eps * (abs(d(m)) + abs(d(m+1)))
                    break
                end
            end
            p = d(l);
            if  m ~= l
                
                if j == 30
                    ierr=l;
                    return
                end
                j = j + 1;
                % form shift
                g = (d(l+1) - p) / (2 * e(l));
                r = sqrt(g*g+1);
                if g  < 0
                    ggg=-abs(r);
                else
                    ggg=abs(r);
                end
                g = d(m) - p + e(l) / (g + ggg);
                s = 1;
                c = 1;
                p = 0;
                mml = m - l;
                % for i=m-1 step -1 until l do --
                for ii = 1: mml
                    i = m - ii;
                    f = s * e(i);
                    b = c * e(i);
                    if abs(f) >= abs(g)
                        c = g / f;
                        r = sqrt(c*c+1.0);
                        e(i+1) = f * r;
                        s = 1.0 / r;
                        c = c * s;
                    else
                        s = f / g;
                        r = sqrt(s*s+1.0);
                        e(i+1) = g * r;
                        c = 1.0 / r;
                        s = s * c;
                    end % if
                    g = d(i+1) - p;
                    r = (d(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d(i+1) = g + p;
                    g = c * r - b;
                    % form first component of vector 
                    f = z(i+1);
                    z(i+1) = s * z(i) + c * f;
                    z(i) = c * z(i) - s * f;
                end
                d(l) = d(l) - p;
                e(l) = g;
                e(m) = 0.0;
            else
                ivar =2;
            end
        end
    end
    % order eigenvalues and eigenvectors 
    for  ii = 2 : n
        i = ii - 1;
        k = i;
        p = d(i);
        for  j = ii : n
            if d(j) <= p
                k = j;
                p = d(j);
            end
        end
        if k ~= i
            d(k) = d(i);
            d(i) = p;
            p = z(i);
            z(i) = z(k);
            z(k) = p;
        end
    end
    % set error -- no convergence to an eigenvalue after 30 iterations 
end
% last card of gausq2 
function[a,b,muzero] = class(kind, n, alpha, beta)
disp('class'),pause,
nargout,pause,
% this procedure supplies the coefficients a(j), b(j) of the recurrence relation
%  b p (x) = (x - a ) p   (x) - b   p   (x)
%  j j            j   j-1       j-1 j-2
% for the various classical (normalized) orthogonal polynomials, and the
% zero-th moment
% muzero = integral w(x) dx of the given polynomial's weight function w(x).
% Since the polynomials are orthonormalized, the tridiagonal matrix is
% guaranteed to be symmetric. The input parameter alpha is used only for
% Laguerre and Jacobi polynomials, and the parameter beta is used only for
% Jacobi polynomials. The Laguerre and Jacobi polynomials require the gamma
% function.
nm1 = n - 1 ;
a = zeros(1,nm1);
b = zeros(1,nm1);
if kind ==1
    % kind = 1: Legendre polynomials p(x) on (-1, +1), w(x) = 1.
    muzero = 2.0;
    for   i = 1: nm1
        a(i) = 0.0;
        abi = i;
        b(i) = abi/sqrt(4*abi*abi - 1.0);
    end
    a(n) = 0.0;
elseif kind==2
    % kind = 2:  Chebyshev polynomials of the first kind t(x)
    % on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
    muzero = pi;
    for   i = 1: nm1
        a(i) = 0.0;
        b(i) = 0.5;
    end
    b(1) = sqrt(0.5) ;
    a(n) = 0.0;
elseif kind==3
    % kind = 3: Chebyshev polynomials of the second kind u(x) on (-1, +1),
    % w(x) = sqrt(1 - x*x)
    muzero = pi/2.0;
    for  i = 1 : nm1
        a(i) = 0.0d0 ;
        b(i) = 0.5d0;
    end
    a(n) = 0.0d0;
elseif kind==4
    % kind = 4: Hermite polynomials h(x) on (-infinity, +infinity),
    % w(x) = exp(-x**2)
    muzero = sqrt(pi);
    for i = 1 : nm1
        a(i) = 0.0;
        b(i) = sqrt(i/2);
    end
    a(n) = 0.0;
elseif kind==5
    % kind = 5: Jacobi polynomials p(alpha, beta)(x) on (-1, +1),
    % w(x) = (1-x)**alpha + (1+x)**beta, alpha and beta greater than -1
    ab = alpha + beta;
    abi = 2.0 + ab;
    muzero = 2.0^(ab+1.)*gamma(alpha+1)*gamma(beta+1) / gamma(abi);
    a(1) = (beta - alpha)/abi;
    b(1) = sqrt(4.0*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi));
    a2b2 = beta*beta - alpha*alpha;
    for  i = 2: nm1
        abi = 2.0*i + ab;
        a(i) = a2b2/((abi - 2.0)*abi);
        b(i) = sqrt(4*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi));
    end
    abi = 2.0*n + ab;
    a(n) = a2b2/((abi - 2.0)*abi);
else
    % kind = 6: Laguerre polynomials l(alpha)(x) on (0, +infinity),
    % w(x) = exp(-x) * x**alpha, alpha greater than -1
    muzero = gamma(alpha + 1.0);
    for   i = 1 :nm1
        a(i) = 2.0*i - 1.0 + alpha;
        b(i) =  sqrt(i*(i + alpha));
    end
    a(n) = 2.0*n - 1 + alpha;
end

function[a,b,muzero] = class1(kind, n, alpha, beta)
% This procedure supplies the coefficients a(j), b(j) of the recurrence
% relation
% b p (x) = (x - a ) p   (x) - b   p   (x)
% j j            j   j-1       j-1 j-2
% for the various classical (normalized) orthogonal polynomials, and the
% zero-th moment muzero = integral w(x) dx of the given polynomial's
% weight function w(x). Since the polynomials are orthonormalized,
% the tridiagonal matrix is guaranteed to be symmetric. The input parameter
% alpha is used only for laguerre and Jacobi polynomials, and the parameter
% beta is used only for Jacobi polynomials. The laguerre and Jacobi
% polynomials require the gamma function.
nm1 = n - 1 ;
a = zeros(1,nm1);
b = zeros(1,nm1);
if kind ==1
    % kind = 1: Legendre polynomials p(x) on (-1, +1), w(x) = 1
    muzero = 2.0;
    for   i = 1: nm1
        a(i) = 0.0;
        abi = i;
        b(i) = abi/sqrt(4*abi*abi - 1.0);
    end
    a(n) = 0.0;
elseif kind==2
    % kind = 2: Chebyshev polynomials of the first kind t(x) on (-1, +1),
    % w(x) = 1 / sqrt(1 - x*x)
    muzero = pi;
    for   i = 1: nm1
        a(i) = 0.0;
        b(i) = 0.5;
    end
    b(1) = sqrt(0.5) ;
    a(n) = 0.0;
elseif kind==3
    % kind = 3: Chebyshev polynomials of the second kind u(x) on (-1, +1),
    % w(x) = sqrt(1 - x*x)
    muzero = pi/2.0;
    for  i = 1 : nm1
        a(i) = 0.0d0 ;
        b(i) = 0.5d0;
    end
    a(n) = 0.0d0;
elseif kind==4
    % kind = 4: Hermite polynomials h(x) on (-infinity, +infinity),
    % w(x) = exp(-x**2)
    muzero = sqrt(pi);
    for i = 1 : nm1
        a(i) = 0.0;
        b(i) = sqrt(i/2);
    end
    a(n) = 0.0;
elseif kind==5
    % kind = 5: Jacobi polynomials p(alpha, beta)(x) on (-1, +1),
    % w(x) = (1-x)**alpha + (1+x)**beta, alpha and beta greater than -1
    ab = alpha + beta;
    abi = 2.0 + ab;
    muzero = 2.0^(ab+1.)*gamma(alpha+1)*gamma(beta+1) / gamma(abi);
    a(1) = (beta - alpha)/abi;
    b(1) = sqrt(4.0*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi));
    a2b2 = beta*beta - alpha*alpha;
    for  i = 2: nm1
        abi = 2.0*i + ab;
        a(i) = a2b2/((abi - 2.0)*abi);
        b(i) = sqrt(4*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi));
    end
    abi = 2.0*n + ab;
    a(n) = a2b2/((abi - 2.0)*abi);
else
    % kind = 6: Laguerre polynomials l(alpha)(x) on (0, +infinity),
    % w(x) = exp(-x) * x**alpha, alpha greater than -1
    muzero = gamma(alpha + 1.0);
    for   i = 1 :nm1
        a(i) = 2.0*i - 1.0 + alpha;
        b(i) =  sqrt(i*(i + alpha));
    end
    a(n) = 2.0*n - 1 + alpha;
end

function s=solve(shift, n, a, b)
% This procedure performs elimination to solve for the n-th component of
% the solution delta to the equation (jn - shift*identity) * delta  = en,
% where en is the vector of all zeroes except for 1 in the n-th position.
% The matrix jn is symmetric tridiagonal, with diagonal elements a(i),
% off-diagonal elements b(i). This equation must be solved to obtain the
% appropriate changes in the lower 2 by 2 submatrix of coefficients for
% orthogonal polynomials.
alpha = a(1) - shift;
nm1 = n - 1;
for i = 2: nm1
    alpha = a(i) - shift - b(i-1)^2/alpha;
    s = 1.0d0/alpha;
end

