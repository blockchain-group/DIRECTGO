function y = hs035(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs035.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 0.111111115111110
%   x* = (1.33333332055107, 0.777777777397707, 0.444444442025613) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+x(2)+2*x(3) <= 3;
%         0 <= x(1) <= 3;
%         0 <= x(2) <= 3;
%         0 <= x(3) <= 3;
%   
% Problem Properties:
%   n  = 3;
%   #g = 1;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 9-8*x(1)-6*x(2)-4*x(3)+2*x(1)^2+2*x(2)^2+x(3)^2+2*x(1)*x(2)+2*x(1)*x(3);
end