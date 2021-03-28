function [c, ceq] = s224c( x )
c(1) = -x(1)-3*x(2);
c(2) = -18+x(1)+3*x(2);
c(3) = -x(1)-x(2); 
c(4) = -8+x(1)+x(2);
ceq = [];
end

