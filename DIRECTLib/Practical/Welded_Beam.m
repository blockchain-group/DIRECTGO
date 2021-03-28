function y = Welded_Beam(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Welded_Beam.m
%
% Original source: 
% - L. C. Cagnina, S. C. Esquivel and C. A. Coello Coello: Solving 
%   Engineering Optimization Problems with the Simple Constrained Particle 
%   Swarm Optimizer, Informatica (Slovenia), Vol.3, No.32, pp. 319-326, 2008. 
%
% Globally optimal solution:
%   f* = 1.72488430718736
%   x* = (0.205555555555556,3.46081305356568,9.07164699373966,
%         0.205555555555556) 
%
% Constraints (including variable bounds):
%   g(1): 0.125-x(1)                                     <= 0;
%   g(2): t-t_c                                          <= 0;
%   g(3): d-d_max                                        <= 0;
%   g(4): 0.10471*x(1)^2+0.04811*x(3)*x(4)*(14+x(2))-5   <= 0;
%   g(5): x(1)-x(4)                                      <= 0;
%   g(6): s-s_max                                        <= 0;
%   g(7): t-t_max                                        <= 0;
%   where
%   P = 6000; L = 14; E = 30e+6; G = 12e+6;
%   t_max = 13600; s_max = 30000; d_max = 0.25;
%   P = 6000; L = 14; E = 30e+6; G = 12e+6;
%   M = P*(L+(x(2)/2)); 
%   J = 2*(sqrt(2)*x(1)*x(2)*(x(2)^2/12+0.25*(x(1)+x(3))^2));
%   P_c = (4.013*E/(6*L^2))*x(3)*x(4)^3*(1-0.25*x(3)*sqrt(E/G)/L);
%   R = sqrt(((x(2)^2)/4)+((x(1)+x(3))/2)^2);
%   t1 = P/(sqrt(2)*x(1)*x(2)); t2 = M*R/J;
%   t = sqrt(t1^2+t1*t2*x(2)/R+t2^2);
%   s = 6*P*L/(x(4)*x(3)^2);
%   d = 4*P*L^3/(E*x(4)*x(3)^3);
%         0.1 <= x(1) <= 2;
%         0.1 <= x(2) <= 10;
%         0.1 <= x(3) <= 10;
%         0.1 <= x(4) <= 2;
%   
% Problem Properties:
%   n  = 4;
%   #g = 7;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 1.10471*x(1)^2*x(2)+0.04811*x(3)*x(4)*(14.0+x(2));
end