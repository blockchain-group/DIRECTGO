function y = Genocop10(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Genocop10.m
%
% Original source: 
% - Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures=
%   Evolution Programs' third edition 1996, Appendix C case 11 pp.223-240.
%
% Globally optimal solution:
%   f* = -4.52836567655105
%   x* = (0.00316227034547004, 2, 5, 1) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(3) <= 4;
%   g(2): x(2)+2*x(4) <= 4;
%         0 <= x(1) <= 3;
%         0 <= x(2) <= 10;
%         0 <= x(3) <= 10;
%         0 <= x(4) <= 1;
%   
% Problem Properties:
%   n  = 4;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -(x(1)^0.6+x(2)^0.6-6*x(1)-4*x(3)+3*x(4));
end