function y = G10(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G010.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = 7049.24802052867
%   x* = (579.306685017979589, 1359.97067807935605, 5109.97065743133317,
%         182.01769963061534,  295.601173702746792, 217.982300369384632,
%         286.41652592786852,  395.601173702746735)
%
% Constraints (including variable bounds):
%   g(1): -1+0.0025*(x(4)+x(6))                       <= 0;
%   g(2): -1+0.0025*(-x(4)+x(5)+x(7))                 <= 0;
%   g(3): -1+0.01*(-x(5)+x(8))                        <= 0;
%   g(4): 100*x(1)-x(1)*x(6)+833.33252*x(4)-83333.333 <= 0;
%   g(5): x(2)*x(4)-x(2)*x(7)-1250*x(4)+1250*x(5)     <= 0;
%   g(6): x(3)*x(5)-x(3)*x(8)-2500*x(5)+1250000       <= 0;
%         100  <= x(1) <= 10000;
%         1000 <= x(2) <= 10000;
%         1000 <= x(3) <= 10000;
%         10   <= x(4) <= 1000;
%         10   <= x(5) <= 1000;
%         10   <= x(6) <= 1000;
%         10   <= x(7) <= 1000;
%         10   <= x(8) <= 1000;
%   
% Problem Properties:
%   n  = 8;
%   #g = 6;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)+x(2)+x(3);  
end