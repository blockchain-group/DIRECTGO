function y = s231(x)
% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------ 
y = 100*((x(2)-x(1)^2)^2)+(1-x(1))^2;
end