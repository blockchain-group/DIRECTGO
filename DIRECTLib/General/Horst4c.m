function [c, ceq] = Horst4c( x )
c(1) = x(1)+x(2)+2*x(3)-6; 
c(2) = x(1)+(1/2)*x(2)-2;
c(3) = -x(2)-2*x(3)+1;
c(4) = -x(1)+(1/2);
ceq = [];
end

