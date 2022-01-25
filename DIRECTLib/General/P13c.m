function [c, ceq] = P13c( x )
ceq(1) = abs(600*x(1) - 50*x(3) - x(1)*x(3) + 5000); 
ceq(2) = abs(600*x(2) + 50*x(3) - 15000); 
c = [];
end
