function y = Shekel7(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Shekel7.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2354.htm
%
% Globally optimal solution:
%   f = -10.4029405668187;
%   x = [4; 4; 4; 4];
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...4
%   bounds = ones(4, 1).*[0, 10];
%   
% Problem Properties:
%   n  = 4;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
a = [4, 1, 8, 6, 3, 2, 5;
     4, 1, 8, 6, 7, 9, 5;
     4, 1, 8, 6, 3, 2, 3;
     4, 1, 8, 6, 7, 9, 3];
c = [0.1 0.2 0.2 0.4 0.4 0.6 0.3];
if size(x, 1) == 1
 x = x';
end
for i = 1:7
 b = (x - a(:, i)).^2;
 d(i) = sum(b);
end
y = -sum((c + d).^(-1));
end