function [c, ceq] = Horst5c( x )
c(1) = x(1)+x(2)+x(3)-2;
c(2) = x(1)+x(2)-(1/4)*x(3)-1;
c(3) = -2*x(1)-2*x(2)+x(3)-1;
c(4) = x(3)-3;
ceq = [];
end

