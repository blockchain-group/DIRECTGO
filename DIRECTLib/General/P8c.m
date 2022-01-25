function [c, ceq] = P8c( x )
c(1) = x(2) - x(1)^2 - 2*x(1) + 2;
c(2) = -x(1) + x(2) - 8;
ceq = [];
end
