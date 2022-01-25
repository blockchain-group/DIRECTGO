function y = P5(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P5.m
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
% Test problem P05 after reformulation contains 2 variables, 2 equality 
% and 2 inequality constraints. In the original problem formulation there 
% were 3 variables and 3 equality constraints.
%
% Globally optimal solution:
%   f* = 201.1593340582
%   x* = (6.29342997676684; 3.82183908126620) 
%
% Constraints (including variable bounds):
%   g(1): (0.5*(x(1)+x(2))^2+150)-267.42                <= 0;
%   g(2): -(0.5*(x(1)+x(2))^2+150)                      <= 0;
%   h(1): 30*x(1)-6*x(1)^2-(0.5*(x(1)+x(2))^2+150)+250   = 0;
%   h(2): 20*x(2)-12*x(2)^2-(0.5*(x(1)+x(2))^2+150)+300  = 0;
%         0 <= x(1) <= 9.422;
%         0 <= x(2) <= 5.903;
%   
% Problem Properties:
%   n  = 2;
%   #g = 2;
%   #h = 2;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 2;
    y.nh = 2;
    y.xl = @(i) 0;
    xu = [9.422, 5.903];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 201.1593340582003009;
    xmin = [6.2934299767668431; 3.8218390812661962];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P5c(i);
    return
end
y = 0.5*((x(1) + x(2))^2) + 150;  
end

function [c, ceq] = P5c( x )
c(1) = (0.5*(x(1) + x(2))^2 + 150) - 267.42; 
c(2) = -(0.5*(x(1) + x(2))^2 + 150); 
ceq(1) = abs(30*x(1) - 6*x(1)^2 - (0.5*(x(1) + x(2))^2 + 150) + 250); 
ceq(2) = abs(20*x(2) - 12*x(2)^2 - (0.5*(x(1) + x(2))^2 + 150) + 300);  
end