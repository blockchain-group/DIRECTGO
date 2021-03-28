function [c, ceq] = G12c( x)
for p=1:9
    for q=1:9
        for r=1:9
            z(p,q,r) = (x(1)-p)^2+(x(2)-q)^2+(x(3)-r)^2-0.0625;
        end
    end
end
for p=1:9
    for q=1:9
        Z1(p,q) = min(z(p,q,:));    
    end
    Z2(p) = min(Z1(p,:));
end
c = min(Z2);

ceq = [];
end