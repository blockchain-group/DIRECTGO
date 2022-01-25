function y = P3a(x)
% -------------------------------------------------------------------------
% MATLAB coding by: Linas Stripinis
% Name:
%   P3a.m
%
% Original source: 
% - Christodoulos A. Floudas, Panos M. Pardalos, Claire S. Adjiman, 
%   William R. Esposito, Zeynep H. Gumus, Stephen T. Harding, 
%   John L. Klepeis, Clifford A. Meyer, Carl A. Schweiger. 1999. Handbook 
%   of Test Problems in Local and Global Optimization. Nonconvex 
%   Optimization and Its Applications, Vol. 33. Springer Science Business 
%   Media, B.V. https://doi.org/10.1007/978-1-4757-3040-1
%
% Globally optimal solution:
%   f* = -0.3888098393593030
%   x* = [0.768423023726840; 0.517221595629136; 0.206569293270163;
%         0.388809839359303; 3.08904144680576; 5.02850614351362] 
%
% Constraints (including variable bounds):
%   g(1): x(5)^(1/2)+x(6)^(1/2)-4                 <= 0;
%   h(1): x(1)+0.09755988*x(1)*x(5)-1              = 0;
%   h(2): x(2)-x(1)+0.0965842812*x(2)*x(6)         = 0;
%   h(3): x(3)+x(1)+0.03919080*x(3)*x(5)-1         = 0;
%   h(4): x(4)-x(3)+x(2)-x(1)+0.03527172*x(4)*x(6) = 0;
%         0       <= x(1) <= 1;
%         0       <= x(2) <= 1;
%         0       <= x(3) <= 1;
%         0       <= x(4) <= 1;
%         10^(-5) <= x(5) <= 16;
%         10^(-5) <= x(6) <= 16;
%   
% Problem Properties:
%   n  = 6;
%   #g = 1;
%   #h = 4;  
% -------------------------------------------------------------------------
if nargin == 0
    y.nx = 6;
    y.ng = 1;
    y.nh = 4;
    xl = [0, 0, 0, 0, 10^(-5), 10^(-5)];
    y.xl = @(i) xl(i);
    xu = [1, 1, 1, 1, 16, 16];
    y.xu = @(i) xu(i);
    y.fmin = @(i) -0.3888098393593030;
    xmin = [0.768423023726840; 0.517221595629136; 0.206569293270163;...
        0.388809839359303; 3.08904144680576; 5.02850614351362];
    y.xmin = @(i) xmin(i);
    y.confun = @(i) P3ac(i);
    return
end
y = -x(4);
end

function [c, ceq] = P3ac(x)
c = x(5)^(1/2) + x(6)^(1/2) - 4; 
ceq(1) = abs(x(1) + 0.09755988*x(1)*x(5)-1);
ceq(2) = abs(x(2) - x(1) + 0.0965842812*x(2)*x(6)); 
ceq(3) = abs(x(3) + x(1) + 0.03919080*x(3)*x(5) - 1); 
ceq(4) = abs(x(4) - x(3) + x(2) - x(1) + 0.03527172*x(4)*x(6));
end