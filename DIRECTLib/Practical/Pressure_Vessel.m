function y = Pressure_Vessel(x)
% -------------------------------------------------------------------------
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
%   f* = 7163.7395688773249275982379913330078125
%   x* = (1.1, 0.625, 56.9948186528497, 51.0012517339584) 
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
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 4;
    y.ng = 6;
    y.nh = 0;
    bounds = [1, 1.375; 0.625, 1; 25, 150; 25, 240];
    y.xl = @(i) bounds(i, 1);
    y.xu = @(i) bounds(i, 2);
    y.fmin = @(i) 7163.7395688773249275982379913330078125;
    xmin = [1.1; 0.625; 56.9948186528497; 51.0012517339584];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Pressure_Vesselc(i);
    return
end
y = 0.6224*x(1)*x(3)*x(4) + 1.7781*x(2)*x(3)^2 +...
    3.1661*x(1)^2*x(4) + 19.84*x(1)^2*x(3); 
end

function [c, ceq] = Pressure_Vesselc( x )
c(1) = -x(1) + 0.0193*x(3);
c(2) = -x(2) + 0.00954*x(3);  
c(3) = -pi*x(3)^2*x(4) - (4/3)*pi*x(3)^3 + 1296000; 
c(4) = x(4) - 240; 
c(5) = 1.1 - x(1);
c(6) = 0.6 - x(2);
ceq=[];
end