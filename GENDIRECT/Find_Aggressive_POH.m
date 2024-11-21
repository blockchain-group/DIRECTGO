function boxes = Find_Aggressive_POH( fc, C, II, Dat, ~ )

% Initialization
Sz = length( C );
m_set = [ fc; Dat.MSS.DD( II ); 1:length( II ) ];

if ( Dat.STA.Eqs )

    % Select all equal candidates
    POH = cell( 1, Sz );
    for i = 1:Sz

        MB = m_set( :, m_set(2, :) == C(i) );
        MV = min( MB(1, :) );
        TM = MB(1, :) == MV;
        POH{i} = MB(3, TM);

    end
    boxes = cell2mat( POH );

else

    % Select only one candidates from several equally
    if ( Dat.STA.Dbc )
        LFmin = Dat.MSS.LFmin(:, II);
    end

    POH = zeros( 3, Sz );
    for i = 1:Sz

        MB = m_set( :, m_set(2, :) == C(i) );
        [ MV, MI ] = min( MB(1, :) );
        TM = find( MB(1, :) == MV );

        if ( length( TM ) ~= 1 )

            if ( ~Dat.STA.Dbc )
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
        POH(:, i) = MB(:, MI);

    end
    POH = flip( POH, 2 );
    boxes = sort( POH(3, :) );

end

return