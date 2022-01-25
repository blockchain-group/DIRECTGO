function y = P2d(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P2d_mod.m
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
% Test problem P2d after reformulation contain 5 variables and 12 
% inequality constraints. In the original problem formulation there were 
% 10 variables, 5 equality and 2 inequality constraints.
%
% Globally optimal solution:
%   f* = -400
%   x* = (0, 100, 0, 100, 1) 
%
% Constraints (including variable bounds):
%   g(1): x(5)*x(1)+2*x(3)-2.5*(x(1) + x(3))                   <= 0;
%   g(2): x(5)*x(2)+2*x(4)-1.5*(x(2)+x(4))                     <= 0;
%   g(3): x(3)+x(4)-300                                        <= 0;
%   g(4): -(x(3)+x(4))                                         <= 0;
%   g(5): x(2)+x(4) - 200                                      <= 0;
%   g(6): -(x(2)+x(4))                                         <= 0;
%   g(7): x(1) + x(3) - 100                                    <= 0;
%   g(8): -(x(1) + x(3))                                       <= 0;
%   g(9): ((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2) -300             <= 0;
%   g(10): -((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2)                <= 0;
%   g(11): (x(1)+x(2)-((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2))-300 <= 0;
%   g(12): -((x(1)+x(2)-((x(1)*x(5)+x(2)*x(5)-x(1)-x(2))/2)))  <= 0;
%         0 <= x(1) <= 100;
%         0 <= x(2) <= 200;
%         0 <= x(3) <= 100;
%         0 <= x(4) <= 200;
%         1 <= x(5) <= 3;
%   
% Problem Properties:
%   n  = 5;
%   #g = 12;
%   #h = 0;  
% ------------------------------------------------------------------------- 
if nargin == 0
    y.nx = 5;
    y.ng = 12;
    y.nh = 0;
    xl = [0, 0, 0, 0, 1];
    y.xl = @(i) xl(i);
    xu = [100, 200, 100, 200, 3];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -400;
    xmin = [0, 100, 0, 100, 1];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P2dc(i);
    return
end
y = -9*(x(1) + x(3)) - 15*(x(2) + x(4)) + 6*((x(5)*(x(1) + x(2)) -...
    x(1) - x(2))/2) + 16*((x(1) + x(2) - ((x(5)*(x(1) + x(2)) -...
    x(1) - x(2))/2))) + 10*(x(3) + x(4)); 
end

function [c, ceq] = P2dc( x )
c(1) = x(5)*x(1) + 2*x(3) - 2.5*(x(1) + x(3));   
c(2) = x(5)*x(2) + 2*x(4) - 1.5*(x(2) + x(4));  
c(3) = x(3) + x(4) - 300; 
c(4) = -(x(3) + x(4));  
c(5) = x(2) + x(4) - 200; 
c(6) = -(x(2) + x(4));  
c(7) = x(1) + x(3) - 100; 
c(8) = -(x(1) + x(3));  
c(9) = ((x(1)*x(5) + x(2)*x(5) - x(1) - x(2))/2) - 300; 
c(10) = -((x(1)*x(5) + x(2)*x(5) - x(1) - x(2))/2);  
c(11) = (x(1) + x(2) - ((x(1)*x(5) + x(2)*x(5) - x(1) - x(2))/2)) - 300;  
c(12) = -((x(1) + x(2) - ((x(1)*x(5) + x(2)*x(5) - x(1) - x(2))/2)));  
ceq = [];
end