function y = zy2(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   zy2.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 2
%   x* = (0; 2; 0) 
%
% Constraints (including variable bounds):
%     g(1) = -x(1)^2 - x(2)^2 - x(3)^2 + 4 <=0  
%     g(2) = x(1)^2 + x(2)^2 + x(3)^2 - 10 <=0  
%     g(3) = x(3) - 5 <=0 
%
%       0 <= x(i) <= 10, i = 1...n
%       bounds = ones(n, 1).*[0, 10];
%   
% Problem Properties:
%   n  = 3;
%   #g = 3;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 2;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 10;
    y.fmin = @(i) 2;
    xmin = [0; 2; 0];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) zy2c(i);
    return
end
    y = x(1)^3 - 6*x(1)^2 + 11*x(1) + x(2) + x(3);
end

function [Ineq, eq] = zy2c(x)
    Ineq(1) = -x(1)^2 - x(2)^2 - x(3)^2 + 4;  
    Ineq(2) = x(1)^2 + x(2)^2 + x(3)^2 - 10; 
    Ineq(3) = x(3) - 5; 
    eq = [];
end