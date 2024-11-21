function V = VertexGeneration(n)

ok = 1;
b = zeros(1, n);
k = 1;
V(:, k) = b';
k = k + 1;
while (ok == 1)
    j = n;
    while ((j > 0) && (b(j) == 1))
        b(j) = 0;
        j = j - 1;
    end
    if (j > 0)
        b(j) = 1;
        V(:, k) = b';
        k = k + 1;
    else
        ok = 0;
    end
end

return