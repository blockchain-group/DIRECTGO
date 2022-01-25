function y = Drop_wave(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Drop_wave.m
%
% Original source:
%  - http://www.sfu.ca/~ssurjano/
%
% Globally optimal solution:
%   f = -1
%   x(i) = [0] i = 1...n;
% 
% Variable bounds:
%   -5.12 <= x(i) <= 6.12, i = 1...n
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 2;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) -5.12;
    y.xu = @(i) +6.12;
    y.fmin = @(i) -1;
    y.xmin = @(i) 0;
    return
end
frac1 = 1 + cos(12*sqrt(x(1)^2 + x(2)^2));
frac2 = 0.5*(x(1)^2 + x(2)^2) + 2;
y = -frac1/frac2;
end