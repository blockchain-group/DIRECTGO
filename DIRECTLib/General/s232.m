function y = s232(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   s232.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -1
%   x* = (3, 1.7320508075688772) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)/sqrt(3) + x(2)     <= 0;
%   g(2): -x(1) - sqrt(3)*x(2)     <= 0;
%   g(3): -6 + x(1) + sqrt(3)*x(2) <= 0;
%         0 <= x(1) <= 100;
%         0 <= x(2) <= 100;
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
    y.xu = @(i) 100;
    y.fmin = @(i) -1;
    xmin = [3, 1.7320508075688772];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) s232c(i);
    return
end
y = -1*(9 - (x(1) - 3)^2)*(x(2)^3/(27*sqrt(3)));
end

function [c, ceq] = s232c( x )
c(1) = -x(1)/sqrt(3) + x(2);
c(2) = -x(1) - sqrt(3)*x(2);
c(3) = -6 + x(1) + sqrt(3)*x(2);
ceq = [];
end