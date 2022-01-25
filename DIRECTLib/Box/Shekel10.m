function y = Shekel10(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Shekel10.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2354.htm
%
% Globally optimal solution:
%   f = -10.5364098166920463
%   x = [4.0007465305280281; 4.0005929353320706; 3.9996634007540983; 
%        3.9995097988662054];
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
    y.fmin = @(i) -10.5364098166920463;
    xmin = [4.0007465305280281; 4.0005929353320706; 3.9996634007540983;... 
            3.9995097988662054];
    y.xmin = @(i) xmin(i);
    return
end
a = [4, 1, 8, 6, 3, 2, 5, 8, 6, 7;
     4, 1, 8, 6, 7, 9, 5, 1, 2, 3.6;
     4, 1, 8, 6, 3, 2, 3, 8, 6, 7;
     4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6];
c = [ 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5];
if size(x, 1) == 1
 x = x';
end
for i = 1:10
 b = (x - a(:, i)).^2;
 d(i) = sum(b);
end
y = -sum((c + d).^(-1));
end