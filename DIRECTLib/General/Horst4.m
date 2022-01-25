function y = Horst4(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Horst4.m
%
% Original source: 
% - Horst, R., Pardalos, P.M., Thoai, N.V. (1995). Introduction to  
%   Global Optimization. Nonconvex Optimization and Its Application. 
%   Kluwer, Dordrecht  
%
% Globally optimal solution:
%   f* = -6.085806194501845
%   x* = (2, 0, 2) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)            <= -(1/2);
%   g(2): -x(2)-2*x(3)     <= -1;
%   g(3): x(1)+(1/2)*x(2)  <= 2;
%   g(4): x(1)+x(2)+2*x(3) <= 6;
%         0 <= x(1) <= 2;
%         0 <= x(2) <= 3;
%         0 <= x(3) <= 2.8;
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
    y.xl = @(i) 0;
    xu = [2, 3, 2.8];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -6.085806194501845;
    xmin = [2, 0, 2];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Horst4c(i);
    return
end
y = -(abs(x(1) + (1/2)*x(2) + (2/3)*x(3)))^(3/2);
end

function [c, ceq] = Horst4c( x )
c(1) = x(1) + x(2) + 2*x(3) - 6; 
c(2) = x(1) + (1/2)*x(2) - 2;
c(3) = -x(2) - 2*x(3) + 1;
c(4) = -x(1) + (1/2);
ceq = [];
end