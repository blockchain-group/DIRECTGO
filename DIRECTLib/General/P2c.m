function y = P2c(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P2c.m
%
% Original source: 
% - Christodoulos A. Floudas, Panos M. Pardalos, Claire S. Adjiman, 
%   William R. Esposito, Zeynep H. Gumus, Stephen T. Harding, 
%   John L. Klepeis, Clifford A. Meyer, Carl A. Schweiger. 1999. Handbook 
%   of Test Problems in Local and Global Optimization. Nonconvex 
%   Optimization and Its Applications, Vol. 33. Springer Science Business 
%   Media, B.V. https://doi.org/10.1007/978-1-4757-3040-1
%
% Problem have been reformulated by some algebraic manipulation aiming to 
% reduce the number of variables and equality constraints.
% - Costa, M. F. P., Rocha, A. M. A. C., & Fernandes, E. M. G. P.  
%   Filter-based DIRECT method for constrained global optimization. 
%   Journal of Global Optimization, 71(3), 517–536. (2018) 
%
% Test problem P2c after reformulation contain 5 variables and 10 
% inequality constraints. In the original problem formulation there were 
% 9 variables, 4 equality and 2 inequality constraints.
%
% Globally optimal solution:
%   f* = -750
%   x* = (0, 0, 1.5, 0, 200) 
%
% Constraints (including variable bounds):
%   g(1): x(4)+x(1)-100                                        <= 0;
%   g(2): -(x(4)+x(1))                                         <= 0;
%   g(3): x(5)+x(2)-200                                        <= 0;
%   g(4): -(x(5)+x(2))                                         <= 0;
%   g(5): x(3)*x(5) + 2*x(2) - 1.5*(x(5)+x(2))                 <= 0;
%   g(6): x(3)*x(4) + 2*x(1) - 2.5*(x(4)+x(1))                 <= 0;
%   g(7): ((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2)-500              <= 0;
%   g(8): -((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2)                 <= 0;
%   g(9): x(4)+x(5) - ((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2) -500 <= 0;
%   g(10): -(x(4)+x(5) - ((x(3)*x(4)+x(3)*x(5)-x(4)-x(5))/2))  <= 0;
%         0 <= x(1) <= 500;
%         0 <= x(2) <= 500;
%         0 <= x(3) <= 500;
%         0 <= x(4) <= 500;
%         0 <= x(5) <= 500;
%   
% Problem Properties:
%   n  = 5;
%   #g = 10;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 5;
    y.ng = 10;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 500;
    y.fmin = @(i) -750;
    xmin = [0, 0, 1.5, 0, 200];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P2cc(i);
    return
end
y = -9*(x(4) + x(1)) - 15*(x(5) + x(2)) + 6*(((x(3)*x(4) +...
    x(3)*x(5) - x(4) - x(5))/2)) + 13*(x(4) + x(5) -...
    ((x(3)*x(4) + x(3)*x(5) - x(4) - x(5))/2)) + 10*(x(1) + x(2)); 
end

function [c, ceq] = P2cc( x )
c(1) = x(4) + x(1) - 600;
c(2) = -(x(4) + x(1)); 
c(3) = x(5) + x(2) - 200; 
c(4) = -(x(5) + x(2));  
c(5) = x(3)*x(5) + 2*x(2) - 1.5*(x(5) + x(2)); 
c(6) = x(3)*x(4) + 2*x(1) - 2.5*(x(4) + x(1));  
c(7) = ((x(3)*x(4) + x(3)*x(5) - x(4) - x(5))/2) - 500; 
c(8) = -((x(3)*x(4) + x(3)*x(5) - x(4) - x(5))/2);  
c(9) = x(4) + x(5) - ((x(3)*x(4) + x(3)*x(5) - x(4) - x(5))/2) - 500; 
c(10) = -(x(4) + x(5) - ((x(3)*x(4) + x(3)*x(5) - x(4) - x(5))/2));  
ceq = [];
end

