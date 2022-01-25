function y = P15(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P15.m
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
%   f* = 0
%   x* = (10.6018948261553; 31.8056843775809; 7.59242078959846) 
%
% Constraints (including variable bounds):
%   h(1): x(1)+x(2)+x(3)-50             = 0;
%   h(2): x(2)/x(1)-3                   = 0;
%   h(3): x(3)^2/(x(1)*x(2)^3)-0.000169 = 0;
%         10^(-5) <= x(1) <= 12.5;
%         10^(-5) <= x(2) <= 37.5;
%         0       <= x(3) <= 50;
%   
% Problem Properties:
%   n  = 3;
%   #g = 0;
%   #h = 3;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 0;
    y.nh = 3;
    xl = [10^(-5), 10^(-5), 0];
    y.xl = @(i) xl(i);
    xu = [12.5, 37.5, 50];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 0;
    xmin = [10.6018948261553; 31.8056843775809; 7.59242078959846];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P15c(i);
    return
end
y = 0; 
end

function [c, ceq] = P15c( x )
c = []; 
ceq(1) = abs(x(3)^2/(x(1)*x(2)^3) - 0.000169); 
ceq(2) = abs(x(2)/x(1) - 3); 
ceq(3) = abs(x(1) + x(2) + x(3) - 50); 
end