function y = P13(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P13.m
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
%   f* = 189.3465728929122349
%   x* = (0.0000100000000010; 16.6666583333016725; 100.0001000001799696) 
%
% Constraints (including variable bounds):
%   h(1): 600*x(1)-50*x(3)-x(1)*x(3)+5000 = 0;
%   h(2): 600*x(2)+50*x(3)-15000          = 0;
%         10^(-5) <= x(1) <= 34;
%         10^(-5) <= x(2) <= 17;
%         100     <= x(3) <= 300;
%   
% Problem Properties:
%   n  = 3;
%   #g = 0;
%   #h = 2;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 0;
    y.nh = 2;
    xl = [10^(-5), 10^(-5), 0];
    y.xl = @(i) xl(i);
    xu = [34, 17, 300];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 189.3465728929122349;
    xmin = [0.0000100000000010; 16.6666583333016725; 100.0001000001799696];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P13c(i);
    return
end
y = 35*x(1)^(0.6) + 35*x(2)^(0.6); 
end

function [c, ceq] = P13c( x )
ceq(1) = abs(600*x(1) - 50*x(3) - x(1)*x(3) + 5000); 
ceq(2) = abs(600*x(2) + 50*x(3) - 15000); 
c = [];
end