function y = Pressure_Vessel(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Pressure_Vessel.m
%
% Original source: 
% - L. C. Cagnina, S. C. Esquivel and C. A. Coello Coello: Solving 
%   Engineering Optimization Problems with the Simple Constrained Particle 
%   Swarm Optimizer, Informatica (Slovenia), Vol.3, No.32, pp. 319-326, 2008. 
%
% Globally optimal solution:
%   f* = 7163.73956887529
%   x* = (1.1, 0.625, 56.994818652849, 51) 
%
% Constraints (including variable bounds):
%   g(1): -x(1)+0.0193*x(3)                       <= 0;
%   g(2): -x(2)+0.00954*x(3)                      <= 0;
%   g(3): -pi*x(3)^2*x(4)-(4/3)*pi*x(3)^3+1296000 <= 0;
%   g(4): x(4)-240                                <= 0;
%   g(5): 1.1-x(1)                                <= 0;
%   g(6): 0.6-x(2)                                <= 0;
%         1     <= x(1) <= 1.375;
%         0.625 <= x(2) <= 1;
%         25    <= x(3) <= 150;
%         25    <= x(4) <= 240;
%   
% Problem Properties:
%   n  = 4;
%   #g = 6;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 0.6224*x(1)*x(3)*x(4)+1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3); 
end