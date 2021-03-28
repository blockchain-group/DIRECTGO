function [Ineq, eq] = s365modc(x)
    P = sqrt(x(2)^2) + x(3)^2;
    Q = sqrt(x(3)^2) + (x(2) - x(1))^2;
    Ineq(1) = -((x(4) - x(6))^2 + (x(5) - x(7))^2 - 4);  
    Ineq(2) = -((x(3)*x(4) - x(2)*x(5))/P - 1); 
    Ineq(3) = -((x(3)*x(6) - x(2)*x(7))/P - 1);
    Ineq(4) = -((x(1)*x(3) + (x(2) - x(1))*x(5) - x(3)*x(4))/Q - 1);
    Ineq(5) = -((x(1)*x(3) + (x(2) - x(1))*x(7) - x(3)*x(6))/Q - 1);
    Ineq(6) = 0.5 - x(1);
    Ineq(7) = 0.5 - x(3);
    Ineq(8) = 1 - x(5);
    Ineq(9) = 1 - x(7);
    eq = [];
end