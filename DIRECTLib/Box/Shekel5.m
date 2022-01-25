function y = Shekel5(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Shekel5.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2354.htm
%
% Globally optimal solution:
%   f = -10.1531996790582
%   x = [4.0000371516773017; 4.0001332773882963; 4.0000371526332925;
%        4.0001332766447479];
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
    y.fmin = @(i) -10.1531996790582308;
    xmin = [4.0000371516773017; 4.0001332773882963;...
            4.0000371526332925; 4.0001332766447479];
    y.xmin = @(i) xmin(i);
    return
end
a = [4, 1, 8, 6, 3;
     4, 1, 8, 6, 7;
     4, 1, 8, 6, 3;
     4, 1, 8, 6, 7];
c = [0.1, 0.2, 0.2, 0.4, 0.4];
if size(x, 1) == 1
 x = x';
end
for i = 1:5
 b = (x - a(:, i)).^2;
 d(i) = sum(b);
end
y = -sum((c + d).^(-1));
end