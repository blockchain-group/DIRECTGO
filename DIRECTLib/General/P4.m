function y = P4(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P4.m
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
%   f* = -20/3
%   x* = (6, 2/3) 
%
% Constraints (including variable bounds):
%   g(1): x(1)*x(2)-4 <= 0;
%         0 <= x(1) <= 6;
%         0 <= x(2) <= 4;
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
    xu = [6, 4];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -20/3;
    xmin = [6, 2/3];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P4c(i);
    return
end
y = -x(1) - x(2); 
end

function [c, ceq] = P4c( x )
c   = x(1)*x(2) - 4; 
ceq = [];
end