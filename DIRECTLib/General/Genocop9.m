function y = Genocop9(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Genocop9.m
%
% Original source: 
% - Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures=
%   Evolution Programs' third edition 1996, Appendix C case 11 pp.223-240.
%
% Globally optimal solution:
%   f* = -2.47142857142857153
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
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 5;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 3;
    y.fmin = @(i) -2.47142857142857153;
    xmin = [1, 0, 0];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Genocop9c(i);
    return
end
y = -((3*x(1) + x(2) - 2*x(3) + 0.8)/(2*x(1) - x(2) + x(3)) + (4*x(1) -...
    2*x(2) + x(3))/(7*x(1) + 3*x(2) - x(3)));
end

function [c, ceq] = Genocop9c(x)
c(1) = x(1) + x(2) - x(3) - 1; 
c(2) = -x(1) + x(2) - x(3) + 1; 
c(3) = 12*x(1) + 5*x(2) + 12*x(2) - 34.8;
c(4) = 12*x(1) + 12*x(2) + 7*x(3) - 29.1;
c(5) = -6*x(1) + x(2) + x(3) + 4.1; 
ceq = [];
end