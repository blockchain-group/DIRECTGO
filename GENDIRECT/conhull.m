function h = conhull(x, y)

    x = x(:); y = y(:);
    xyAR = [x y];
    xyUnique = unique(xyAR, 'rows');
    x = xyUnique(:, 1);
    y = xyUnique(:, 2);
    m = length(x);
    Z = cell(1, m);
    xy = [x, y];
    if m == 2
        for i = 1:m
            Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
        end
        h = cell2mat(Z);
        return;
    end
    if m == 1
        for i = 1:m
            Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
        end
        h = cell2mat(Z);
        return;
    end
    START = 1;
    v = START;
    w = length(x);
    flag = 0;
    h = (1:length(x))';
    while (mod(v, m) + 1 ~= START) || (flag==0)
        flag = mod(v, m) + 1 == w;
        a = v;
        b = mod(v, m) + 1;
        c = mod(mod(v, m) + 1, m) + 1;
        leftturn = sign(det([x(a) y(a) 1; x(b) y(b) 1; x(c) y(c) 1])) >= 0;
        if leftturn
            v = mod(v, m) + 1;
        else
            j = mod(v, m) + 1;
            x(j) = [];
            y(j) = [];
            h(j) = [];
            m = m - 1;
            w = w - 1;
            v = mod(v - 2, m) + 1;
        end
    end
    xy = [x, y];
    hh = size(xy, 1);
    Z = cell(1, hh);
    for i = 1:hh
        Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
    end
    h = cell2mat(Z);
return