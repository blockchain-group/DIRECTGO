function y = G12(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G12.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = -1.0
%   x* = (5,5,5)
%
% Constraints (including variable bounds):
%   g(1): (x(1)-P)^2+(x(2)-Q)^2+(x(3)-R)^2-0.0625 <= 0; where P=[1:9]; R=[1:9]; Q=[1:9]; 
%         0 <= x(1) <= 100;
%         0 <= x(2) <= 100;
%         0 <= x(3) <= 100;
%   
% Problem Properties:
%   n  = 3;
%   #g = 1;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 1;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 100;
    y.fmin = @(i) -1;
    xmin = [5, 5, 5];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) G12c(i);
    return
end
y = -(100 - (x(1) - 5)^2 - (x(2) - 5)^2 - (x(3) - 5)^2)/100;
end

function [c, ceq] = G12c( x)
for p = 1:9
    for q = 1:9
        for r = 1:9
            z(p, q, r) = (x(1) - p)^2 + (x(2)-q)^2 + (x(3) - r)^2 - 0.0625;
        end
    end
end
for p = 1:9
    for q = 1:9
        Z1(p, q) = min(z(p, q, :));    
    end
    Z2(p) = min(Z1(p, :));
end
c = min(Z2);
ceq = [];
end