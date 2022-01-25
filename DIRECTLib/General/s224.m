function y = s224(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   s224.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -304
%   x* = (4, 4)
%
% Constraints (including variable bounds):
%   g(1): -x(1)-3*x(2)    <= 0;
%   g(2): -18+x(1)+3*x(2) <= 0;
%   g(3): -x(1)-x(2)      <= 0;
%   g(4): -8+x(1)+x(2)    <= 0;
%         0 <= x(1) <= 6;
%         0 <= x(2) <= 6;
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
    y.xl = @(i) 0;
    y.xu = @(i) 6;
    y.fmin = @(i) -304;
    xmin = [4, 4];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) s224c(i);
    return
end
y = 2*(x(1)^2) + (x(2)^2) - 48*x(1) - 40*x(2);
end

function [c, ceq] = s224c( x )
c(1) = -x(1) - 3*x(2);
c(2) = -18 + x(1) + 3*x(2);
c(3) = -x(1) - x(2); 
c(4) = -8 + x(1) + x(2);
ceq = [];
end