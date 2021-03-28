function Value = s365mod(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   s365mod.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = 52.139904360199466
%   x* = [1], i = 1...n
%
% Constraints (including variable bounds):
%     P = sqrt(x(2)^2) + x(3)^2;
%     Q = sqrt(x(3)^2) + (x(2) - x(1))^2;
%     g(1) = -((x(4) - x(6))^2 + (x(5) - x(7))^2 - 4) <=0  
%     g(2) = -((x(3)*x(4) - x(2)*x(5))/P - 1) <=0  
%     g(3) = -((x(3)*x(6) - x(2)*x(7))/P - 1) <=0 
%     g(4) = -((x(1)*x(3) + (x(2) - x(1))*x(5) - x(3)*x(4))/Q - 1) <=0 
%     g(5) = -((x(1)*x(3) + (x(2) - x(1))*x(7) - x(3)*x(6))/Q - 1) <=0 
%     g(6) = 0.5 - x(1) <=0 
%     g(7) = 0.5 - x(3) <=0 
%     g(8) = 1 - x(5) <=0 
%     g(9) = 1 - x(7) <=0 
%
%       0 <= x(i) <= 19, i = 1...n
%       bounds = ones(n, 1).*[0, 19];
%   
% Problem Properties:
%   n  = 7;
%   #g = 9;
%   #h = 0;  
% ------------------------------------------------------------------------------ 
    Value = x(1)*x(3);
end

