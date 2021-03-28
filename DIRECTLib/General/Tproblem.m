function f = Tproblem(x)
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
%   x* = (-) 
%
% Constraints (including variable bounds):
%   g(1): sum(x(i)^2) - n <= 0;, i = 1...n
%         -i <= x(i) <= i;, i = 1...n
%   
% Problem Properties:
%   n  = any dimension;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
n = length(x);
f = 0;
for i = 1:n
    f = f + x(i);
end
end