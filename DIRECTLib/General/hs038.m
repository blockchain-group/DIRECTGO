function y = hs038(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs038.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 0
%   x* = (1, 1, 1, 1) 
%
% Constraints (including variable bounds):
%   g(1): x(1)+2*x(2)+2*x(3)  <= 72;
%   g(2): -x(1)-2*x(2)-2*x(3) <= 0;
%         -10 <= x(1) <= 10;
%         -10 <= x(2) <= 10;
%         -10 <= x(3) <= 10;
%         -10 <= x(4) <= 10;
%   
% Problem Properties:
%   n  = 4;
%   #g = 2;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
y = 100*((x(2)-x(1)^2)^2)+((1-x(1))^2)+90*((x(4)-x(3)^2)^2)+((1-x(3))^2)+10.1*(((x(2)-1)^2)+((x(4)-1)^2))+19.8*(x(2)-1)*(x(4)-1);
end