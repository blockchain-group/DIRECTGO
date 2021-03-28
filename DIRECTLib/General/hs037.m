function y = hs037(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs037.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -3456
%   x* = (24, 12, 12) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(2)+2*x(3)  <= 72;
%   g(2): -x(1)-2*x(2)-2*x(3) <= 0;
%         0 <= x(1) <= 42;
%         0 <= x(2) <= 42;
%         0 <= x(3) <= 42;
%   
% Problem Properties:
%   n  = 3;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)*x(2)*x(3);
end