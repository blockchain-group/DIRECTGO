function y = hs076(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs076.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -4.6818181818181818
%   x* = (3/11, 23/11,0, 6/11) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(2)+x(3)+x(4)   <= 5;
%   g(2): 3*x(1)+x(2)+2*x(3)-x(4) <= 4;
%   g(3): x(2)+4*x(3)             >= 1.5;
%         0 <= x(1) <= 1;
%         0 <= x(2) <= 3;
%         0 <= x(3) <= 1;
%         0 <= x(4) <= 1;
%   
% Problem Properties:
%   n  = 4;
%   #g = 3;
%   #h = 0;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 4;
    y.ng = 3;
    y.nh = 0;
    y.xl = @(i) 0;
    xu = [1, 3, 1, 1];
    y.xu = @(i) 42;
    y.fmin = @(i) -4.6818181818181818;
    xmin = [3/11, 23/11,0, 6/11];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) hs076c(i);
    return
end
y = (x(1)^2) + 0.5*(x(2)^2) + (x(3)^2) + 0.5*x(4)^2 -...
    x(1)*x(3) + x(3)*x(4) - x(1) - 3*x(2) + x(3) - x(4);
end

function [c, ceq] = hs076c( x )
c(1) = x(1) + 2*x(2) + x(3) + x(4) - 5;
c(2) = 3*x(1) + x(2) + 2*x(3) - x(4) - 4;
c(3) = -x(2) - 4*x(3) + 1.5;
ceq = [];
end