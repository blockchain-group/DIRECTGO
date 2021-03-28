function y = Hartman3(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Hartman3.m
%
% Original source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1488.htm
%
% Globally optimal solution:
%   f = -3.86278214780000
%   x = [0.114614; 0.555649; 0.852547]
%
% Variable bounds:
%   0 <= x(i) <= 1, i = 1...3
%   bounds = ones(3, 1).*[0, 1];
%   
% Problem Properties:
%   n  = 3;
%   #g = 0;
%   #h = 0;
% ------------------------------------------------------------------------------
a(:, 2) = 10.0*ones(4, 1);
for j=1:2
    a(2*j - 1, 1) = 3.0;
    a(2*j, 1) = 0.1; 
    a(2*j - 1, 3) = 30.0; 
    a(2*j, 3) = 35.0;
end
c(1) = 1.0; 
c(2) = 1.2; 
c(3) = 3.0; 
c(4) = 3.2;
p(1, 1) = 0.36890; 
p(1, 2) = 0.11700; 
p(1, 3) = 0.26730;
p(2, 1) = 0.46990; 
p(2, 2) = 0.43870; 
p(2, 3) = 0.74700;
p(3, 1) = 0.10910; 
p(3, 2) = 0.87320; 
p(3, 3) = 0.55470;
p(4, 1) = 0.03815; 
p(4, 2) = 0.57430; 
p(4, 3) = 0.88280;
s = 0;
for i = 1:4
    sm = 0;
    for j = 1:3
        sm = sm + a(i, j)*(x(j) - p(i, j))^2;
    end
    s = s + c(i)*exp(-sm);
end
y = -s;
end