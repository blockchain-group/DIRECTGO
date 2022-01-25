function y = Tension_Compression_Spring(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Tension_Compression_Spring.m
%
% Original source: 
% - L. C. Cagnina, S. C. Esquivel and C. A. Coello Coello: Solving 
%   Engineering Optimization Problems with the Simple Constrained Particle 
%   Swarm Optimizer, Informatica (Slovenia), Vol.3, No.32, pp. 319-326, 2008. 
%
% Globally optimal solution:
%   f* = 0.01267867559268556071350
%   x* = (0.05169590656, 0.35688327343, 11.2933789329) 
%
% Constraints (including variable bounds):
%   g(1): 1-(((x(2)^3)*x(3))/(71875*x(1)^4))                                          <= 0;
%   g(2): (((4*x(2)-x(1))*x(2))/(12566*(x(1)^3)*(x(2)-x(1))))+(2.46/(12566*x(1)^2))-1 <= 0;
%   g(3): 1-((140.54*x(1))/(x(3)*x(2)^2))                                             <= 0;
%   g(4): ((x(1)+x(2))/1.5)-1                                                         <= 0;
%         0.05 <= x(1) <= 0.2;
%         0.25 <= x(2) <= 1.3;
%         2    <= x(3) <= 15;
%   
% Problem Properties:
%   n  = 3;
%   #g = 4;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 4;
    y.nh = 0;
    bounds = [0.05, 0.2; 0.25, 1.3; 2, 15];
    y.xl = @(i) bounds(i, 1);
    y.xu = @(i) bounds(i, 2);
    y.fmin = @(i) 0.01267867559268556071350;
    xmin = [0.05169590656; 0.35688327343; 11.2933789329];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Tension_Compression_Springc(i);
    return
end
y = x(1)^2*x(2)*(x(3) + 2);
end

function [c, ceq] = Tension_Compression_Springc( x )
c(1) = 1 - (((x(2)^3)*x(3))/(71875*x(1)^4));
c(2) = (((4*x(2)-x(1))*x(2))/(12566*(x(1)^3)*(x(2) - x(1)))) +...
    (2.46/(12566*x(1)^2)) - 1;
c(3) = 1 - ((140.54*x(1))/(x(3)*x(2)^2));
c(4) = ((x(1) + x(2))/1.5) - 1; 
ceq = [];
end