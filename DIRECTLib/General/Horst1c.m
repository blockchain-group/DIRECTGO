function [c, ceq] = Horst1c( x )
c(1) = -4*x(1)+2*x(2)-1; 
c(2) = x(1)+x(2)-4;
c(3) = x(1)-4*x(2)-1;
ceq = [];
end

