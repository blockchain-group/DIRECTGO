function y = hs024(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs024.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -1
%   x* = (3, 1.732050807568876971131999) 
%
% Constraints (including variable bounds):
%   g(1): x(1)/sqrt(3)-x(2)  >= 0;
%   g(2): x(1)+sqrt(3)*x(2)  >= 0;
%   g(3): -x(1)-sqrt(3)*x(2) >= -6;
%         0 <= x(1) <= 5;
%         0 <= x(2) <= 5;
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
    y.xu = @(i) 5;
    y.fmin = @(i) -1;
    xmin = [3, 1.732050807568876971131999];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) hs024c(i);
    return
end
y = (((x(1) - 3)^2) - 9)*((x(2)^3)/(27*sqrt(3)));
end

function [c, ceq] = hs024c( x )
c(1) = -x(1)/sqrt(3) + x(2);
c(2) = -x(1) - sqrt(3)*x(2); 
c(3) = x(1) + sqrt(3)*x(2) - 6; 
ceq = [];
end