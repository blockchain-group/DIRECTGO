function S = Find_RedPareto_POH( fc, C, II, Dat, i_min )

% Initialization
szes = Dat.MSS.DD(II);
if ( Dat.STA.Dbc )
    LFmin = Dat.MSS.LFmin(:, II);
end
sz = length(C);
d_min = cell2mat(arrayfun(@(i) min(fc(szes == C(i))), 1:sz, 'UniformOutput', false));

Su = cell(1, sz);
if Dat.STA.Eqs
    Su = arrayfun(@(i) find((fc == d_min(i)) & (szes == C(i))), 1:sz, 'UniformOutput', false);
else
    if ( ~Dat.STA.Dbc )
        Su = arrayfun(@(i) find((fc == d_min(i)) & (szes == C(i)), 1, 'first'), 1:sz, 'UniformOutput', false);
    else
        for i = 1:sz
            TM = find((fc == d_min(i)) & (szes == C(i)));
            if size(TM, 2) ~= 1
                ii = 0;
                while length(LFmin(:, 1)) ~= ii
                    ii  = ii + 1;
                    if ~isnan(min(LFmin(ii, TM))) && length(TM) > 1
                        TM = TM(LFmin(ii, TM) == min(LFmin(ii, TM)));
                    else
                        break;
                    end
                end
                if length(TM) > 1
                    ii = find(~isnan(LFmin(:, TM(1))), 1, 'last');
                    TM = TM(find(LFmin(ii, TM) == min(LFmin(ii, TM)), 1, 'last'));
                end
                [~, MI] = max(TM);
                TM = TM(MI);
            end
            Su{i} = TM;
        end
    end
end
S_1 = cell2mat( Su );
Su = cell( 1, sz );
d_length = length( C );
if ( d_length > 0 )
    a1 = szes(i_min);
    b1 = fc(i_min);
    a2 = C(d_length);
    b2 = d_min(d_length);
    slope = (b2 - b1) / (a2 - a1);
    const = b1 - slope * a1;
    for i = 1 : length( S_1 )
        j = S_1(i);
        if ( fc(j) <= slope * szes(j) + const + 1E-6 )
            Su{i} = j;
        end
    end
    S_2 = cell2mat( Su );
    if ( isempty( S_2 ) )
        S_2 = S_1;
    end
    xx = szes(S_2);
    yy = fc(S_2);
    h = conhull( xx, yy );
    S_3 = S_2(h);
else
    S_3 = S_1;
end

if length( S_3 ) > 2
    i1 = szes(S_3) == szes(S_3(1));
    i2 = szes(S_3) == szes(S_3(end));
    S = [S_3(i1), S_3(i2)];
else
    S = S_3;
end

return