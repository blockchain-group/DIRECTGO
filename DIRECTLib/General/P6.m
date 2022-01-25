function y = P6(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P6.m
%
% Original source: 
% Original source: 
% - Christodoulos A. Floudas, Panos M. Pardalos, Claire S. Adjiman, 
%   William R. Esposito, Zeynep H. Gumus, Stephen T. Harding, 
%   John L. Klepeis, Clifford A. Meyer, Carl A. Schweiger. 1999. Handbook 
%   of Test Problems in Local and Global Optimization. Nonconvex 
%   Optimization and Its Applications, Vol. 33. Springer Science Business 
%   Media, B.V. https://doi.org/10.1007/978-1-4757-3040-1
%
% Globally optimal solution:
%   f* = 376.2919323266
%   x* = [8.1700182298244304; 7.5607442427640450]
%
% Constraints (including variable bounds):
%   g(1): -x(1) + ((0.2458*x(1)^2)/x(2)) + 6 <= 0;
%         0       <= x(1) <= 115.8;
%         10^(-5) <= x(2) <= 30;
%   
% Problem Properties:
%   n  = 2;
%   #g = 1;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 1;
    y.nh = 0;
    xl = [0, 10^(-5)];
    y.xl = @(i) xl(i);
    xu = [115.8, 30];
    y.xu = @(i) xu(i);
    y.fmin = @(i) 376.2919323265910521;
    xmin = [8.1700182298244304; 7.5607442427640450];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P6c(i);
    return
end
y = 29.4*x(1) + 18*x(2); 
end

function [c, ceq] = P6c( x )
c   = -x(1) + ((0.2458*x(1)^2)/x(2)) + 6; 
ceq = [];
end