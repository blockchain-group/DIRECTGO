function y = P11(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P11.m
%
% Original source: 
% - Christodoulos A. Floudas, Panos M. Pardalos, Claire S. Adjiman, 
%   William R. Esposito, Zeynep H. Gumus, Stephen T. Harding, 
%   John L. Klepeis, Clifford A. Meyer, Carl A. Schweiger. 1999. Handbook 
%   of Test Problems in Local and Global Optimization. Nonconvex 
%   Optimization and Its Applications, Vol. 33. Springer Science Business 
%   Media, B.V. https://doi.org/10.1007/978-1-4757-3040-1
%
% Globally optimal solution:
%   f* = -0.5
%   x* = (0.5, 0.5) 
%
% Constraints (including variable bounds):
%   g(1): 4*x(1)*x(2)+2*x(1)+2*x(2)-3 <= 0;
%         0 <= x(1) <= 2;
%         0 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 1;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 2;
    y.fmin = @(i) -0.5;
    xmin = [0.5, 0.5];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P11c(i);
    return
end
y = -2*x(1)*x(2); 
end

function [c, ceq] = P11c( x )
c   = 4*x(1)*x(2) + 2*x(1) + 2*x(2) - 3; 
ceq = [];
end