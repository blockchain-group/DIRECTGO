function [c, ceq] = s250c( x )
c(1) = -x(1)-2*x(2)-2*x(3); 
c(2) = -72+x(1)+2*x(2)+2*x(3);
ceq = [];
end

