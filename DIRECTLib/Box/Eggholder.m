function y = Eggholder(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Eggholder.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/egg.html
%
% Globally optimal solution:
%   f = -959.6406627208516511
%   x = [512.0000000000002273; 404.2318050882936404]
%
% Variable bounds:
%   -512 <= x(i) <= 512, i = 1...2
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
    y.xl = @(i) -512;
    y.xu = @(i) +512;
    y.fmin = @(i) -959.6406627208516511;
    xmin = [512.0000000000002273; 404.2318050882936404];
    y.xmin = @(i) xmin(i);
    return
end
term1 = -(x(2)+47)*sin(sqrt(abs(x(2) + x(1)/2 + 47)));
term2 = -x(1)*sin(sqrt(abs(x(1) - (x(2) + 47))));
y = term1 + term2;
end