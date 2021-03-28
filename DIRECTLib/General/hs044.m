function y = hs044(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs044.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -15
%   x* = (0, 3, 0, 4) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(2)   <= 8;
%   g(2): 4*x(1)+x(2)   <= 12;
%   g(3): 3*x(1)+4*x(2) <= 12;
%   g(4): 2*x(3)+x(4)   <= 8;
%   g(5): x(3)+2*x(4)   <= 8;
%   g(6): x(3)+x(4)<=5;
%         0 <= x(1) <= 42;
%         0 <= x(2) <= 42;
%         0 <= x(3) <= 42;
%         0 <= x(4) <= 42;
%   
% Problem Properties:
%   n  = 4;
%   #g = 6;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = x(1)-x(2)-x(3)-x(1)*x(3)+x(1)*x(4)+x(2)*x(3)-x(2)*x(4);
end