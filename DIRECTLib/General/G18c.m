function [c, ceq] = G18c( x )
c(1) = x(3)^2 + x(4)^2 - 1; 
c(2) = x(9)^2 - 1;
c(3) = x(5)^2 + x(6)^2 - 1;  
c(4) = x(1)^2 + (x(2) - x(9))^2 - 1; 
c(5) = (x(1) - x(5))^2 + (x(2) - x(6))^2 - 1; 
c(6) = (x(1) - x(7))^2 + (x(2) - x(8))^2 - 1;
c(7) = (x(3) - x(5))^2 + (x(4) - x(6))^2 - 1; 
c(8) = (x(3) - x(7))^2 + (x(4) - x(8))^2 - 1; 
c(9) = (x(7))^2 + (x(8) - x(9))^2 - 1; 
c(10) = x(2)*x(3) - x(1)*x(4);
c(11) = - x(3)*x(9); 
c(12) = x(5)*x(9);
c(13) = x(6)*x(7) - x(5)*x(8); 
ceq=[];
end
