function [c, ceq] = s232c( x )
c(1) = -x(1)/sqrt(3)+x(2);
c(2) = -x(1)-sqrt(3)*x(2);
c(3) = -6+x(1)+sqrt(3)*x(2);
ceq = [];
end

