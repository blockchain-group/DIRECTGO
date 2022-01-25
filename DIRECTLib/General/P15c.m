function [c, ceq] = P15c( x )
c = []; 
ceq(1) = abs(x(3)^2/(x(1)*x(2)^3) - 0.000169); 
ceq(2) = abs(x(2)/x(1) - 3); 
ceq(3) = abs(x(1) + x(2) + x(3) - 50); 
end
