function y = Non_Reg_n3_T100(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   Non_Reg_n3_T100.m
%
% Original source: 
% - Gillard, J.W., Kvasov, D.E., 2017. Lipschitz optimization methods for
%   fitting a sum of damped sinusoids to a series of observations. Statistics
%   and its Interface doi:10.4310/SII.2017.v10.n1.a6.. 
%
% Globally optimal solution:
%   f* = 0
%   x* = [-0.2; 0.4; 0.3]
%
% Box constraints:
%   -1 <= x(1) <= 0;
%   0  <= x(2) <= 1;
%   0  <= x(3) <= 1;
%   
% Problem Properties:
%   n  = 3;
%   #g = 0;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 3;
    y.ng = 0;
    y.nh = 0;
    bounds = [-1, 0; 0, 1; 0, 1];
    y.xl = @(i) bounds(i, 1);
    y.xu = @(i) bounds(i, 2);
    y.fmin = @(i) 0;
    xmin = [-0.2; 0.4; 0.3];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) Three_bar_trussc(i);
    return
end
T = 100; 
zhig_y4 = zeros(1,T);
for j=1:T
  zhig_y4(j) = 1*exp(-0.2*j)*sin(2*pi*0.4*j + 0.3);
end
xt = zeros(1,T);
for j=1:T
   xt(j) = exp(x(1)*j)*sin(2*pi*x(2)*j + x(3));
end
y = 0;
for j=1:T
  zz_temp = (zhig_y4(j) - xt(j))^2;
  y = y + zz_temp;
end