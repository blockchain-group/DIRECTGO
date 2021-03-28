function y = Cross_in_Tray(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Cross_in_Tray.m
%
% Original source:
%  - https://www.sfu.ca/~ssurjano/crossit.html
%
% Globally optimal solution:
%   f = -2.0626118708 
%   x(i) = [1.3491], i = 1...2
%   x = ones(2, 1)*1.3491;
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...2
%   bounds = ones(2, 1).*[0, 10];
%
% Problem Properties:
%   n  = 2;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
fact1 = sin(x(1))*sin(x(2));
fact2 = exp(abs(100 - sqrt(x(1)^2 + x(2)^2)/pi));
y = -0.0001*(abs(fact1*fact2) + 1)^0.1;
end