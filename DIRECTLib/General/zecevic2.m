function y = zecevic2(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   zecevic2.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -4.1250
%   x* = (1.75, 0.25) 
%
% Constraints (including variable bounds):
%     g(1) = x(1) + x(2) - 2 <=0  
%     g(2) = x(1) + 4*x(2) - 4 <=0  
%
%       0 <= x(i) <= 10, i = 1...n
%       bounds = ones(n, 1).*[0, 10];
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
    y.xl = @(i) 0;
    y.xu = @(i) 10;
    y.fmin = @(i) -4.1250;
    xmin = [1.75, 0.25];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) zecevic2c(i);
    return
end
    y = 2*x(2)^2 - 2*x(1) - 3*x(2);
end

function [Ineq, eq] = zecevic2c(x)
    Ineq(1) = x(1) + x(2) - 2;  
    Ineq(2) = x(1) + 4*x(2) - 4; 
    eq = [];
end