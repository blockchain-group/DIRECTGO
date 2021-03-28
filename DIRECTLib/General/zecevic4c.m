function [Ineq, eq] = zecevic4c( x )
    Ineq(1) = x(1)*x(2) - x(1) - x(2);  
    Ineq(2) = -x(1) - x(2) + 3; 
    eq=[];
end

