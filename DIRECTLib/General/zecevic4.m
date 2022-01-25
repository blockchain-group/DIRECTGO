function y = zecevic4(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   zecevic4.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 7.5575077689317425
%   x* = (4.97095269841094; 1.25182835834412) 
%
% Constraints (including variable bounds):
%     g(1) = x(1)*x(2) - x(1) - x(2) <=0  
%     g(2) = -x(1) - x(2) + 3 <=0  
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
    y.fmin = @(i) 7.5575077689317425;
    xmin = [4.9709529494020677; 1.2518287203958374];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) zecevic4c(i);
    return
end
    y = 6*x(1)^2 + x(2)^2 - 60*x(1) - 8*x(2) + 166;
end

function [Ineq, eq] = zecevic4c( x )
    Ineq(1) = x(1)*x(2) - x(1) - x(2);  
    Ineq(2) = -x(1) - x(2) + 3; 
    eq=[];
end