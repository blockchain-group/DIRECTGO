function y = Genocop11(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Genocop11.m
%
% Original source: 
% - Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures=
%   Evolution Programs' third edition 1996, Appendix C case 11 pp.223-240.
%
% Globally optimal solution:
%   f* = -11
%   x* = (0, 6, 0, 1, 1, 0) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(2)+8*x(3)+x(4)+3*x(5)+5*x(6)          <=16;
%   g(2): -8*x(1)-4*x(2)-2*x(3)+2*x(4)+4*x(5)-x(6)       <=-1;
%   g(3): 2*x(1)+0.5*x(2)+0.2*x(3)-3*x(4)-x(5)-4*x(6)    <=24;
%   g(4): 0.2*x(1)+2*x(2)+0.1*x(3)-4*x(4)+2*x(5)+2*x(6)  <=12;
%   g(5): -0.1*x(1)-0.5*x(2)+2*x(3)+5*x(4)-5*x(5)+3*x(6) <=2; 
%         0 <= x(1) <= 10;
%         0 <= x(2) <= 10;
%         0 <= x(3) <= 10;
%         0 <= x(4) <= 1;
%         0 <= x(5) <= 1;
%         0 <= x(6) <= 2;
%   
% Problem Properties:
%   n  = 6;
%   #g = 5;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 6.5*x(1)-0.5*x(1)^2-x(2)-2*x(3)-3*x(4)-2*x(5)-x(6);
end