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
%   f = -959.6406627209
%   x = [512; 404.2319]
%
% Variable bounds:
%   -512 <= x(i) <= 512, i = 1...2
%   bounds = ones(2, 1).*[-512, 512]
%   
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
term1 = -(x(2)+47)*sin(sqrt(abs(x(2) + x(1)/2 + 47)));
term2 = -x(1)*sin(sqrt(abs(x(1) - (x(2) + 47))));
y = term1 + term2;
end