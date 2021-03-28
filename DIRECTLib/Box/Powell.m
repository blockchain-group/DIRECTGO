function y = Powell(xx)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Powell.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2720.htm
%
% Globally optimal solution:
%   f = 0
%   x(i) = [0], i = 1...n
%   x = zeros(n, 1);
%
% Variable bounds:
%   -4 <= x(i) <= 5, i = 1...n
%   bounds = ones(n, 1).*[-4, 5];
%   
% Problem Properties:
%   n  = 4;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
d = length(xx);
sum = 0;

for ii = 1:(d/4)
	term1 = (xx(4*ii-3) + 10*xx(4*ii-2))^2;
	term2 = 5 * (xx(4*ii-1) - xx(4*ii))^2;
	term3 = (xx(4*ii-2) - 2*xx(4*ii-1))^4;
	term4 = 10 * (xx(4*ii-3) - xx(4*ii))^4;
	sum = sum + term1 + term2 + term3 + term4;
end

y = sum;
end