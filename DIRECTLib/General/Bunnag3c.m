function [c, ceq] = Bunnag3c( x )
c(1) = x(1) + 2*x(4) - 4; 
c(2) = 3*x(1) + 3*x(4) + x(5) - 4; 
c(3) = 2*x(2) + 4*x(4) + 2*x(5) - 6; 
ceq = [];
end

