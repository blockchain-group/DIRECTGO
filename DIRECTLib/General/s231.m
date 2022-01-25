function y = s231(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   s231.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 0
%   x* = (1, 1)
%
% Constraints (including variable bounds):
%   g(1): x(1)/3 - x(2) - 0.1  <= 0;
%   g(2): -x(1)/3 - x(2) - 0.1 <= 0;
%         -10 <= x(1) <= 10;
%         -10 <= x(2) <= 10;
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
    y.xl = @(i) -10;
    y.xu = @(i) 10;
    y.fmin = @(i) 0;
    xmin = [1, 1];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) s231c(i);
    return
end
y = 100*((x(2) - x(1)^2)^2) + (1 - x(1))^2;
end

function [c, ceq] = s231c( x )
c(1) = -x(1)/3 - x(2) - 0.1;
c(2) = x(1)/3 - x(2) - 0.1;
ceq = [];
end