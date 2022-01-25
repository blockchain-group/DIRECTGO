function y = zecevic3(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   zecevic3.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 97.3094501419674316
%   x* = (2.7955451883416216; 1.0885435682323457) 
%
% Constraints (including variable bounds):
%     g(1) = -x(1)*x(2) + 1 <=0  
%     g(2) = x(1)^2+x(2)^2 - 9 <=0  
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
    y.fmin = @(i) 97.3094501419674316;
    xmin = [2.7955451883416216; 1.0885435682323457];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) zecevic3c(i);
    return
end
    y = 7*x(1)^2 + 3*x(2)^2 - 84*x(1) - 24*x(2) + 300;
end

function [Ineq, eq] = zecevic3c(x)
    Ineq(1) = -x(1)*x(2) + 1;  
    Ineq(2) = x(1)^2+x(2)^2 - 9; 
    eq = [];
end