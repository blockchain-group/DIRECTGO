function [c, ceq] = Tproblemc( x )
n = length(x);
ff = 0;
for i = 1:n
    ff = ff + x(i)^2;
end
c = ff - n;
ceq=[];
end

