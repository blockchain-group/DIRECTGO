function y = G24(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G24.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = -5.50801327159536
%   x* = (2.32952019747762, 3.17849307411774) 
%
% Constraints (including variable bounds):
%   g(1): -2*x(1)^4+8*x(1)^3-8*x(1)^2+x(2)-2            <= 0;
%   g(2): -4*x(1)^4+32*x(1)^3-88*x(1)^2+96*x(1)+x(2)-36 <= 0;
%          0 <= x(1) <= 3;
%          0 <= x(2) <= 4;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)-x(2);
end