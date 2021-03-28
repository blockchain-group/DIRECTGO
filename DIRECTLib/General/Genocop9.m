function y = Genocop9(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Genocop9.m
%
% Original source: 
% - Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures=
%   Evolution Programs' third edition 1996, Appendix C case 11 pp.223-240.
%
% Globally optimal solution:
%   f* = -2.4714
%   x* = (1, 0, 0) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+x(2)-x(3)         <= 1;
%   g(2): -x(1)+x(2)-x(3)        <= -1;
%   g(3): 12*x(1)+5*x(2)+12*x(2) <= 34.8; 
%   g(4): 12*x(1)+12*x(2)+7*x(3) <= 29.1;
%   g(5): -6*x(1)+x(2)+x(3)      <= -4.1;
%         0 <= x(1) <= 3;
%         0 <= x(2) <= 3;
%         0 <= x(3) <= 3;
%   
% Problem Properties:
%   n  = 3;
%   #g = 5;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -((3*x(1)+x(2)-2*x(3)+0.8)/(2*x(1)-x(2)+x(3))+(4*x(1)-2*x(2)+x(3))/(7*x(1)+3*x(2)-x(3)));
end