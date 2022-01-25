function y = P10(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P10.m
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
%   f* = 0.7417819582470551
%   x* = (0.1294095225512604; 0.4829629131445343) 
%
% Constraints (including variable bounds):
%   g(1): -16*x(1)*x(2)+1      <= 0;
%   g(2): -4*x(1)^2-4*x(2)^2+1 <= 0;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 2;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 1;
    y.fmin = @(i) 0.7417819582470551;
    xmin = [0.1294095225512604; 0.4829629131445343];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P10c(i);
    return
end
y = 2*x(1) + x(2); 
end

function [c, ceq] = P10c( x )
c(1) = -16*x(1)*x(2) + 1;
c(2) = -4*x(1)^2 - 4*x(2)^2 + 1;
ceq = [];
end