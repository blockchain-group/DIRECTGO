function y = Horst2(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst2.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -6.899519052838329
%   x* = (2.5, 0.75) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)+x(2)  <= 1;
%   g(2): x(1)-2*x(2) <= 1;
%   g(3): x(1)+2*x(2) <= 4;
%         0 <= x(1) <= 2.5;
%         0 <= x(2) <= 2;
%   
% Problem Properties:
%   n  = 2;
%   #g = 3;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 3;
    y.nh = 0;
    y.xl = @(i) 0;
    xu = [2.5, 2];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -6.899519052838329;
    xmin = [2.5, 0.75];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Horst2c(i);
    return
end
y = -x(1)^2 - x(2)^(3/2);
end

function [c,ceq] = Horst2c( x )
c(1) = x(1)+2*x(2) - 4;
c(2) = x(1)-2*x(2) - 1;
c(3) = -x(1) + x(2) - 1;
ceq = [];
end