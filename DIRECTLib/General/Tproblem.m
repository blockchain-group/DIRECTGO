function y = Tproblem(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Tproblem.m
%
% Original source: 
% - Finkel D. E. (2005). Global Optimization with the DIRECT 
%   Algorithm. Raleigh, North Carolina State University.
%
% Globally optimal solution:
%   f* = -n
%   x* = [-1, -1, ..., -1]  
%
% Constraints (including variable bounds):
%   g(1): sum(x(i)^2) - n <= 0;, i = 1...n
%         -i <= x(i) <= i;, i = 1...n
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 1;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 0;
    y.ng = 1;
    y.nh = 0;
    y.xl = @(i) -i;
    y.xu = @(i) +i;
    y.fmin = @(i) -i;
    y.xmin = @(i) -1;
    y.confun = @(i) Tproblemc(i);
    return
end

n = length(x);
y = 0;
for i = 1:n
    y = y + x(i);
end
end

function [c, ceq] = Tproblemc( x )
n = length(x);
ff = 0;
for i = 1:n
    ff = ff + x(i)^2;
end
c = ff - n;
ceq = [];
end