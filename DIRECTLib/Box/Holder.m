function y = Holder(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Holder.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/holder.html
%
% Globally optimal solution:
%   f = -19.2085025678867538
%   x = [8.0550234733225885, 9.6645900114093131]
%
% Variable bounds:
%   -10 <= x(i) <= 10, i = 1...2
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
    y.xl = @(i) -10;
    y.xu = @(i) +10;
    y.fmin = @(i) -19.2085025678867538;
    xmin = [8.0550234733225885, 9.6645900114093131];
    y.xmin = @(i) xmin(i);
    return
end
fact1 = sin(x(1))*cos(x(2));
fact2 = exp(abs(1 - sqrt(x(1)^2 + x(2)^2)/pi));
y = -abs(fact1*fact2);
end