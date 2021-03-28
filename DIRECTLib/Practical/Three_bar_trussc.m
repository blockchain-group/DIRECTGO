function [c, ceq] = Three_bar_trussc( x )
c(1) = ((sqrt(2)*x(1)+x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*2-2; 
c(2) = ((x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*2-2; 
c(3) = ((1)/(x(1) + sqrt(2)*x(2)))*2-2; 
ceq=[];
end

