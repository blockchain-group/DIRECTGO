function y = hs021(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs021.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -99.96
%   x* = (2, 0) 
%
% Constraints (including variable bounds):
%   g(1): 10*x(1)-x(2) >= 10;
%         2   <= x(1) <= 50;
%         -50 <= x(2) <= 50;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)^2/100+x(2)^2-100;
end