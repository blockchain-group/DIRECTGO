function y = hs024(x)
% ------------------------------------------------------------------------------
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
%   x* = (3, 1.73205080752269) 
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
% ------------------------------------------------------------------------------ 
y = (((x(1)-3)^2)-9)*((x(2)^3)/(27*sqrt(3)));
end