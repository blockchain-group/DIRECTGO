function y = McCormick(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   McCormick.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/mccorm.html
%
% Globally optimal solution:
%   f = -1.9132229549
%   x = [0.54719; 1.54719]
%
% Variable bounds:
%   -1.5 <= x(1) <= 4
%   -3 <= x(2) <= 4
%   bounds = [-1.5 ,4; -3, 4];
%
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
term1 = sin(x(1) + x(2));
term2 = (x(1) - x(2))^2;
term3 = -1.5*x(1);
term4 = 2.5*x(2);
y = term1 + term2 + term3 + term4 + 1;
end