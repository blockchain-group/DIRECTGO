function y = G06(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G06.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = -6961.81387558015
%   x* = (14.09500000000000064, 0.8429607892154795668)
%
% Constraints (including variable bounds):
%   g(1): -(x(1)-5)^2-(x(2)-5)^2+100  <= 0;
%   g(2): (x(1)-6)^2+(x(2)-5)^2-82.81 <= 0;
%          13 <= x(1) <= 100;
%          0  <= x(2) <= 100;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = (x(1)-10)^3+(x(2)-20)^3;  
end