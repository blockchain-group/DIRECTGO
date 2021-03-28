function y = hs036(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs036.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 3300
%   x* = (20, 11, 15) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(2)+2*x(3) <= 72;
%         0 <= x(1) <= 20;
%         0 <= x(2) <= 11;
%         0 <= x(3) <= 15;
%   
% Problem Properties:
%   n  = 3;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = -x(1)*x(2)*x(3);
end