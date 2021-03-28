function [Ineq, eq] = zecevic3c(x)
    Ineq(1) = -x(1)*x(2) + 1;  
    Ineq(2) = x(1)^2+x(2)^2 - 9; 
    eq = [];
end

