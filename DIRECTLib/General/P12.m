function y = P12(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P12.m
%
% Original source: 
% - Christodoulos A. Floudas, Panos M. Pardalos, Claire S. Adjiman, 
%   William R. Esposito, Zeynep H. Gumus, Stephen T. Harding, 
%   John L. Klepeis, Clifford A. Meyer, Carl A. Schweiger. 1999. Handbook 
%   of Test Problems in Local and Global Optimization. Nonconvex 
%   Optimization and Its Applications, Vol. 33. Springer Science Business 
%   Media, B.V. https://doi.org/10.1007/978-1-4757-3040-1
%
% Problem have been reformulated by some algebraic manipulation aiming to 
% reduce the number of variables and equality constraints.
% - Costa, M. F. P., Rocha, A. M. A. C., & Fernandes, E. M. G. P.  
%   Filter-based DIRECT method for constrained global optimization. 
%   Journal of Global Optimization, 71(3), 517–536. (2018) 
%
% Test problem P12 after reformulation contains 1 variable and 2 inequality
% constraints. In the original problem formulation there were 2 variables
% and 1 equality constraints.
%
% Globally optimal solution:
%   f* = -16.738893184394637
%   x* = (0.717536188588019) 
%
% Constraints (including variable bounds):
%   g(1): 2-2*x(1)^4-3  <= 0;
%   g(2): -(2-2*x(1)^4) <= 0;
%         0 <= x(1) <= 2;
%   
% Problem Properties:
%   n  = 1;
%   #g = 2;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 1;
    y.ng = 2;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 2;
    y.fmin = @(i) -16.738893184394637;
    y.xmin = @(i) 0.717536188588019;
    y.confun = @(i) P12c(i);
    return
end
y = -12*x(1) + 6*x(1)^4 + 4*x(1)^8 - 10; 
end

function [c, ceq] = P12c( x )
c(1) = 2 - 2*x(1)^4 - 3; 
c(1) = -(2 - 2*x(1)^4);   
ceq = [];
end