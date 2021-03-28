function [Ineq, eq] = zy2c(x)
    Ineq(1) = -x(1)^2 - x(2)^2 - x(3)^2 + 4;  
    Ineq(2) = x(1)^2 + x(2)^2 + x(3)^2 - 10; 
    Ineq(3) = x(3) - 5; 
    eq = [];
end
