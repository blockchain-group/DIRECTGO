function y = s250(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   s250.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -3300
%   x* = (20, 11, 15)
%
% Constraints (including variable bounds):
%   g(1): -x(1)-2*x(2)-2*x(3)    <= 0;
%   g(2): -72+x(1)+2*x(2)+2*x(3) <= 0;
%         0 <= x(1) <= 20;
%         0 <= x(2) <= 11;
%         0 <= x(3) <= 40;
%   
% Problem Properties:
%   n  = 3;
%   #g = 2;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 2;
    y.nh = 0;
    y.xl = @(i) 0;
    xu = [20, 11, 40];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -3300;
    xmin = [20, 11, 15];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) s250c(i);
    return
end
y = -1*x(1)*x(2)*x(3);
end

function [c, ceq] = s250c( x )
c(1) = -x(1) - 2*x(2) - 2*x(3); 
c(2) = -72 + x(1) + 2*x(2) + 2*x(3);
ceq = [];
end