function boxes = Find_Pareto_POH( fc, C, II, Dat, ~ )

% Initialization
if ( Dat.STA.Dbc ) && ( ~Dat.STA.Eqs )
    LFmin = Dat.MSS.LFmin(:, II);
end
SMS = zeros( 3, length( C ) );
m_set = [fc; Dat.MSS.DD( II ); 1:length( II )];

for i = 1:length( C )
    MB = m_set(:, m_set(2, :) == C(i));
    [MV, MI] = min( MB(1, :) );
    TM = find( MB(1, :) == MV );
    if ( length( TM ) ~= 1 )
        if ( ~Dat.STA.Dbc ) || ( Dat.STA.Eqs )
            [ ~, MI ] = max( MB(3, TM) );
            MI = TM(MI);
        else
            ii = 0;
            while ( length( LFmin(:, 1) ) ~= ii )

                ii = ii + 1;
                if ( ~isnan( min( LFmin(ii, MB(3, TM)) ) ) ) && ( length( TM ) > 1 )
                    TM = TM( LFmin(ii, MB(3, TM)) == min( LFmin(ii, MB(3, TM)) ) );
                else
                    break;
                end

            end
            if ( length( TM ) > 1 )

                pp = MB(3, TM);
                ii = find( ~isnan( LFmin(:, pp(1)) ), 1, 'last' );
                TM = TM(find( LFmin(ii, MB(3, TM)) == min( LFmin(ii, MB(3, TM)) ), 1, 'last' ));

            end
            [ ~, MI ] = max( MB(3, TM) );
            MI = TM(MI);

        end
    end
    SMS(:, i) = MB(:, MI);
end
SMS = flip(SMS,2);
fc_min = SMS(1, :);
s_i = 0;
setas = zeros(1, size(fc_min, 2));
index = size(fc_min, 2);

while ( index ~= 0 )
    [ m_m, index ] = min( fc_min(1:index) );
    if ( ~isnan( m_m ) )
        s_i = s_i + 1;
        setas(s_i) = index;
    end
    index = index - 1;
end
setas = setas(1:s_i);

if ( Dat.STA.Ags )

    % Select all equal candidates
    POH = cell( 1, size( C, 2 ) );
    for i = 1:size( setas, 2 )
        MB = m_set(:, m_set(2, :) == C(setas(i)));
        MV = min(MB(1, :));
        TM = MB(1, :) == MV;
        POH{i} = MB(3, TM);
    end
    boxes = cell2mat(POH);

else
    
    % Select only one candidates from several equally
    boxes = SMS(3, setas);

end

return