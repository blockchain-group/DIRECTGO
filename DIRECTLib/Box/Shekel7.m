function y = Shekel7(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Shekel7.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2354.htm
%
% Globally optimal solution:
%   f = -10.4029405668187;
%   x = [4.0005729159315848; 4.0006893648356527; 3.9994897106343918; 
%        3.9996061608131148];
%
% Variable bounds:
%   0 <= x(i) <= 10, i = 1...4
%   
% Problem Properties:
%   n  = 4;
%   #g = 0;
%   #h = 0;
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 4;
    y.ng = 0;
    y.nh = 0;
    y.xl = @(i) 0;
    y.xu = @(i) 10;
    y.fmin = @(i) -10.4029405668186641;
    xmin = [4.0005729159315848; 4.0006893648356527; 3.9994897106343918;...
            3.9996061608131148];
    y.xmin = @(i) xmin(i);
    return
end
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