function y = P7(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P7.m
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
%   f* = -2.8284271047
%   x* = (-1.41421356338664;-1.41421354136081) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)+x(2)-1     <= 0;
%   g(2): x(1)-x(2)-1      <= 0;
%   g(3): -x(1)^2-x(2)^2+1 <= 0;
%   g(4): x(1)^2+x(2)^2-4  <= 0;
%         -2 <= x(1) <= 2;
%         -2 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 4;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 4;
    y.nh = 0;
    y.xl = @(i) -2;
    y.xu = @(i) 2;
    y.fmin = @(i) -2.8284271247459052;
    xmin = [-1.4142141971953714;-1.4142129275505337];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P7c(i);
    return
end
y = x(1) + x(2); 
end

function [c, ceq] = P7c( x )
c(1) = -x(1) + x(2) - 1; 
c(2) = x(1) - x(2) - 1;
c(3) = -x(1)^2 - x(2)^2 + 1;
c(4) = x(1)^2 + x(2)^2 - 4; 
ceq = [];
end

