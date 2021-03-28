function y = Hartman6(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Hartman6.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1488.htm
%
% Globally optimal solution:
%   f = -3.32236801141551
%   x = [0.20169; 0.150011; 0.476874; 0.275332; 0.311652; 0.6573] 
%
% Variable bounds:
%   0 <= x(i) <= 1, i = 1...6
%   bounds = ones(6, 1).*[0, 1];
%   
% Problem Properties:
%   n  = 6;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
alpha = [1.0, 1.2, 3.0, 3.2]';
A = [10, 3, 17, 3.5, 1.7, 8;
     0.05, 10, 17, 0.1, 8, 14;
     3, 3.5, 1.7, 10, 17, 8;
     17, 8, 0.05, 10, 0.1, 14];
P = 10^(-4) * [1312, 1696, 5569, 124, 8283, 5886;
               2329, 4135, 8307, 3736, 1004, 9991;
               2348, 1451, 3522, 2883, 3047, 6650;
               4047, 8828, 8732, 5743, 1091, 381];
outer = 0;
for ii = 1:4
	inner = 0;
	for jj = 1:6
		xj = x(jj);
		Aij = A(ii, jj);
		Pij = P(ii, jj);
		inner = inner + Aij*(xj - Pij)^2;
	end
	new = alpha(ii)*exp(-inner);
	outer = outer + new;
end

y = -(2.58 + outer)/1.94;
end