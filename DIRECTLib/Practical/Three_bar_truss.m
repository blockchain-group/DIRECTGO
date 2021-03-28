function y = Three_bar_truss(x)
% ------------------------------------------------------------------------------
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
%   f* = 263.8958453583
%   x* = (0.78868, 0.40825) 
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
% ------------------------------------------------------------------------------ 
y = (2*sqrt(2)*x(1)+x(2))*100;
end