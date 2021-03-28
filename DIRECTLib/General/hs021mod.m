function y = hs021mod(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs021mod.m
%
% Original source: 
% - Hock, W., & Schittkowski, K. (1980). Test examples for 
%   nonlinear programming codes. Journal of Optimization Theory and 
%   Applications, 30(1), 127–129. https://doi.org/10.1007/BF00934594
%
% Globally optimal solution:
%   f* = 4.04
%   x* = (2, 10, 0, 2, 0, 0, 0)
%
% Constraints (including variable bounds):
%   g(1): 10*x(1)-x(2)-10 <= 0;
%         2   <= x(1) <= 50;
%         -50 <= x(2) <= 50;
%         0   <= x(3) <= 50;
%         2   <= x(4) <= 10;
%         -10 <= x(5) <= 10;
%         -10 <= x(6) <= 0;
%         0   <= x(7) <= 10;
%   
% Problem Properties:
%   n  = 7;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -100+0.01*(x(1)^2+x(3)^2+x(5)^2+x(6)^2)+(x(2)^2+x(4)^2+x(7)^2);
end