function y = Tension_Compression_Spring(x)
% ------------------------------------------------------------------------------
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
%   f* = 0.0126793170000000
%   x* = (0.05169590656, 0.35688327343, 11.2933789329 ) 
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
% ------------------------------------------------------------------------------ 
y = x(1)^2*x(2)*(x(3)+2);
end