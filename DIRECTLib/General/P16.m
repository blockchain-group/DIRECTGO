function y = P16(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P16.m
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
% Test problem P16 after reformulation contains 2 variables and 6 
% inequality constraints. In the original problem formulation there were 5
% variables and 3 equality constraints.

%
% Globally optimal solution:
%   f* = 0.7049249272475995
%   x* = (1.8201759971679992; 2.9560114983146040) 
%
% Constraints (including variable bounds):
%   g(1): ((x(1)-1)/(36-12*x(1)))-1.5834  <= 0;
%   g(2): ((x(2)-x(1))/(32-8*x(2)))-3.625 <= 0;
%   g(3): ((5-x(2))/4)-1                  <= 0;
%   g(4): -((x(1)-1)/(36-12*x(1)))        <= 0;
%   g(5): -((x(2)-x(1))/(32-8*x(2)))      <= 0;
%   g(6): -((5-x(2))/4)                   <= 0;
%         1 <= x(1) <= 3;
%         1 <= x(2) <= 4;
%   
% Problem Properties:
%   n  = 2;
%   #g = 6;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 6;
    y.nh = 0;
    y.xl = @(i) 1;
    xu =[3, 4];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 0.7049249272475995;
    xmin = [1.8201759971679992; 2.9560114983146040];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P16c(i);
    return
end
y = ((x(1) - 1)/(36 - 12*x(1))) + ((x(2) -...
    x(1))/(32 - 8*x(2))) + ((5 - x(2))/4); 
end

function [c, ceq] = P16c( x )
c(1) = ((x(1) - 1)/(36 - 12*x(1))) - 1.5834; 
c(2) = ((x(2) - x(1))/(32 - 8*x(2))) - 3.625; 
c(3) = ((5 - x(2))/4) - 1; 
c(4) = -((x(1) - 1)/(36 - 12*x(1))); 
c(5) = -((x(2) - x(1))/(32 - 8*x(2))); 
c(6) = -((5 - x(2))/4); 
ceq = [];
end