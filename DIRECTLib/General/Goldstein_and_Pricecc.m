function [c, ceq] = Goldstein_and_Pricecc( x )
c(1) = -((x(1)-1)^2)-((x(2)-1)^2)+0.9;
c(2) = -((x(1)+1)^2)-((x(2)+1)^2)+1.1;
ceq = [];
end

