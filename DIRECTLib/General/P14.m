function y = P14(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P14.m
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
% Test problem P14 after reformulation contains 3 variables and 4 
% inequality constraints. In the original problem formulation there were 4
% variables, 1 equality and 2 inequality constraints.
%
% Globally optimal solution:
%   f* = -4.5142016513619279
%   x* = (4/3, 4, 0) 
%
% Constraints (including variable bounds):
%   g(1): (1/3)*x(2)-x(1)-2          <= 0;
%   g(2): x(1)+2*((1/3)*x(2)-x(1))-4 <= 0;
%   g(3): x(2)+2*x(3)-4              <= 0;
%   g(4): -((1/3)*x(2)-x(1))         <= 0;
%         10^(-5) <= x(1) <= 3;
%         10^(-5) <= x(2) <= 4;
%         0       <= x(3) <= 1;
%   
% Problem Properties:
%   n  = 3;
%   #g = 4;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 4;
    y.nh = 0;
    xl = [10^(-5), 10^(-5), 0];
    y.xl = @(i) xl(i);
    xu = [3, 4, 1];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -4.5142016513619279;
    xmin = [4/3, 4, 0];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P14c(i);
    return
end
y = x(1)^(0.6) + x(2)^(0.6) - 2*x(1) - (4/3)*x(2) + 3*x(3); 
end

function [c, ceq] = P14c( x )
c(1) = (1/3)*x(2) - x(1) - 2; 
c(2) = x(1) + 2*((1/3)*x(2) - x(1)) - 4; 
c(3) = x(2) + 2*x(3) - 4; 
c(4) = -((1/3)*x(2) - x(1)); 
ceq = [];
end