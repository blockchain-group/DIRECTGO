function y = Three_bar_truss(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Three_bar_truss.m
%
% Original source: 
% - L. C. Cagnina, S. C. Esquivel and C. A. Coello Coello: Solving 
%   Engineering Optimization Problems with the Simple Constrained Particle 
%   Swarm Optimizer, Informatica (Slovenia), Vol.3, No.32, pp. 319-326, 2008. 
%
% Globally optimal solution:
%   f* = 263.8958433764917686
%   x* = (0.7886753129194131; 0.4082477860859604) 
%
% Constraints (including variable bounds):
%   g(1): ((sqrt(2)*x(1)+x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*2-2 <= 0;
%   g(2): ((x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*2-2              <= 0;
%   g(3): ((1)/(x(1) + sqrt(2)*x(2)))*2-2                        <= 0;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1;
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
    bounds = [0, 1; 0, 1];
    y.xl = @(i) bounds(i, 1);
    y.xu = @(i) bounds(i, 2);
    y.fmin = @(i) 263.8958433764917686;
    xmin = [0.7886753129194131; 0.4082477860859604];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Three_bar_trussc(i);
    return
end
y = (2*sqrt(2)*x(1) + x(2))*100;
end

function [c, ceq] = Three_bar_trussc( x )
c(1) = ((sqrt(2)*x(1) + x(2))/(sqrt(2)*x(1)^2 + 2*x(1)*x(2)))*2 - 2; 
c(2) = ((x(2))/(sqrt(2)*x(1)^2 + 2*x(1)*x(2)))*2 - 2; 
c(3) = ((1)/(x(1) + sqrt(2)*x(2)))*2 - 2; 
ceq = [];
end