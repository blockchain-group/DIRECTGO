function y = G08(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G08.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = -0.095825041418035855622648
%   x* = (1.22797135260752599, 4.24537336612274885) 
%
% Constraints (including variable bounds):
%   g(1): 1-x(1)+(x(2)-4)^2 <= 0;
%   g(2): x(1)^2-x(2)+1     <= 0;
%          0 <= x(1) <= 100;
%          0 <= x(2) <= 100;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 2;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 100;
    y.fmin = @(i) -0.095825041418035855622648;
    xmin = [1.22797135260752599, 4.24537336612274885];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) G08c(i);
    return
end
y = -(sin(2*pi*x(1))^3*sin(2*pi*x(2)))/(x(1)^3*(x(1) + x(2)));  
end

function [c, ceq] = G08c( x )
c(1) = 1 - x(1) + (x(2) - 4)^2; 
c(2) = x(1)^2 - x(2) + 1;  
ceq=[];
end