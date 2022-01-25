function y = P1(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P1.m
%
% Original source: 
% - Christodoulos A. Floudas, Panos M. Pardalos, Claire S. Adjiman, 
%   William R. Esposito, Zeynep H. Gumus, Stephen T. Harding, 
%   John L. Klepeis, Clifford A. Meyer, Carl A. Schweiger. 1999. Handbook 
%   of Test Problems in Local and Global Optimization. Nonconvex 
%   Optimization and Its Applications, Vol. 33. Springer Science Business 
%   Media, B.V. https://doi.org/10.1007/978-1-4757-3040-1
%
% Globally optimal solution:
%   f* = 0.0293325909344364817166184167263
%   x* = [1.11921973023266, 1.22305228147849, 1.53652051685888,...
%        1.96627013264726, 1.78704855345026]
%
% Constraints (including variable bounds):
%   g(1): x(1)+x(2)^2+x(3)^3-3*sqrt(2)-2 = 0;
%   g(2): x(2)-x(3)^2+x(4)-2*sqrt(2)+2   = 0;
%   g(3): x(1)*x(5)-2                    = 0;
%         -5 <= x(1) <= 5;
%         -5 <= x(2) <= 5;
%         -5 <= x(3) <= 5;
%         -5 <= x(4) <= 5;
%         -5 <= x(5) <= 5;
%   
% Problem Properties:
%   n  = 5;
%   #g = 0;
%   #h = 3;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 5;
    y.ng = 0;
    y.nh = 3;
    y.xl = @(i) -5;
    y.xu = @(i) 5;
    y.fmin = @(i) 0.0293325909344364817166184167263;
    xmin = [1.11921973023266, 1.22305228147849, 1.53652051685888,...
        1.96627013264726, 1.78704855345026];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P1c(i);
    return
end
y = (x(1) - 1)^2 + (x(1) - x(2))^2 + (x(2) - x(3))^3 +...
    (x(3) - x(4))^4 + (x(4) - x(5))^4;
end

function [c, ceq] = P1c( x )
ceq(1) = abs(x(2) - x(3)^2 + x(4) - 2*sqrt(2) + 2); 
ceq(2) = abs(x(1) + x(2)^2 + x(3)^3 - 3*sqrt(2) - 2); 
ceq(3) = abs(x(1)*x(5) - 2) - 10^(-4); 
c = [];
end