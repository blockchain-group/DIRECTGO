function y = hs076(x)
% ------------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   hs076.m
%
% Original source: 
% - Vaz, A.I.F.: PSwarm solver home page (2010).
%   http://www.norg.uminho.pt/aivaz/pswarm/. Accessed 12 Dec 2013  
%
% Globally optimal solution:
%   f* = -4.68181818021818
%   x* = (0.272727275053595,2.09090907978328,4.63157892977425e-10,
%         0.545454563156689) 
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
% ------------------------------------------------------------------------------ 
y = (x(1)^2)+0.5*(x(2)^2)+(x(3)^2)+0.5*x(4)^2-x(1)*x(3)+x(3)*x(4)-x(1)-3*x(2)+x(3)-x(4);
end