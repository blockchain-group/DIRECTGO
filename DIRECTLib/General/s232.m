function y = s232(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   s232.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -0.999999999600000
%   x* = (3, 1.73205080733794) 
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
% ------------------------------------------------------------------------------ 
y = -1*(9-(x(1)-3)^2)*(x(2)^3/(27*sqrt(3)));
end