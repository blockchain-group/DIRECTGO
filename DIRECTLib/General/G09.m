function y = G09(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   G09.m
%
% Original source: 
% - Suganthan, P. N., Hansen, N., Liang, J. J., Deb, K., Chen, Y.-P., 
%   Auger, A., & Tiwari, S. (2005). Problem Definitions and Evaluation 
%   Criteria for the CEC 2006 Special Session on Constrained Real-Parameter
%   Optimization. KanGAL, (May), 251–256. https://doi.org/c
%
% Globally optimal solution:
%   f* = 680.630057374402
%   x* = (2.33049935147405174, 1.95137236847114592,   -0.477541399510615805,
%         4.36572624923625874, -0.624486959100388983, 1.03813099410962173,
%         1.5942266780671519)
%
% Constraints (including variable bounds):
%   g(1): -127+2*x(1)^2+3*x(2)^4+x(3)+4*x(4)^2+5*x(5)         <= 0;
%   g(2): -282+7*x(1)+3*x(2)+10*x(3)^2+x(4)-x(5)              <= 0;
%   g(3): -196+23*x(1)+x(2)^2+6*x(6)^2-8*x(7)                 <= 0;
%   g(4): 4*x(1)^2+x(2)^2-3*x(1)*x(2)+2*x(3)^2+5*x(6)-11*x(7) <= 0;
%         -10 <= x(1) <= 10;
%         -10 <= x(2) <= 10;
%         -10 <= x(3) <= 10;
%         -10 <= x(4) <= 10;
%         -10 <= x(5) <= 10;
%         -10 <= x(6) <= 10;
%         -10 <= x(7) <= 10;
%   
% Problem Properties:
%   n  = 7;
%   #g = 4;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = (x(1)-10)^2+5*(x(2)-12)^2+x(3)^4+3*(x(4)-11)^2+10*x(5)^6+7*x(6)^2+x(7)^4-4*x(6)*x(7)-10*x(6)-8*x(7);  
end