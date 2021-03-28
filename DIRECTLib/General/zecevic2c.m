function [Ineq, eq] = zecevic2c(x)
    Ineq(1) = x(1) + x(2) - 2;  
    Ineq(2) = x(1) + 4*x(2) - 4; 
    eq = [];
end

