function y = G01(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G01.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = -15
%   x* = (1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1)
%
% Constraints (including variable bounds):
%   g(1): 2*x(1)+2*x(2)+x(10)+x(11)-10 <= 0;
%   g(2): 2*x(1)+2*x(3)+x(10)+x(12)-10 <= 0;
%   g(3): 2*x(2)+2*x(3)+x(11)+x(12)-10 <= 0;
%   g(4): -8*x(1)+x(10)                <= 0;
%   g(5): -8*x(2)+x(11)                <= 0;
%   g(6): -8*x(3)+x(12)                <= 0;
%   g(7): -2*x(4)-x(5)+x(10)           <= 0;
%   g(8): -2*x(6)-x(7)+x(11)           <= 0;
%   g(9): -2*x(8)-x(9)+x(12)           <= 0;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 1;
%         0 <= x(3) <= 1;
%         0 <= x(4) <= 1;
%         0 <= x(5) <= 1;
%         0 <= x(6) <= 1;
%         0 <= x(7) <= 1;
%         0 <= x(8) <= 1;
%         0 <= x(9) <= 1;
%         0 <= x(10) <= 100;
%         0 <= x(11) <= 100;
%         0 <= x(12) <= 100;
%         0 <= x(13) <= 1;
%   
% Problem Properties:
%   n  = 13;
%   #g = 9;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 5*sum(x(1:4))-5*sum(x(1:4).^2)-sum(x(5:13));  
end