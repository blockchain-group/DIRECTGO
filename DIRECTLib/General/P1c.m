function [c, ceq] = P1c( x )
ceq(1) = abs(x(2) - x(3)^2 + x(4) - 2*sqrt(2) + 2); 
ceq(2) = abs(x(1) + x(2)^2 + x(3)^3 - 3*sqrt(2) - 2); 
ceq(3) = abs(x(1)*x(5) - 2) - 10^(-4); 
c = [];
end

