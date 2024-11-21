function [ MSS, VAL ] = DBVS( VAL, MSS, POH )

for i = 1:size( POH, 2 )

    Dim = VAL.condition * VAL.n + ~VAL.condition;

    V = MSS.LL(1:VAL.n, MSS.SS(POH(i)).V);
    longDistD = roundn( pdist2( V', V' ), -16 );
    [LG, kk] = max( max( longDistD, [], 1 ) );
    [~, ll] = max( longDistD(:, kk) );

    for t = 1:Dim

        % New simplex
        VAL.I = VAL.I + 1;

        % Partitioning point on the longest side
        midPoint = (V(:, kk) + V(:, ll)) / 2;

        % Partitioning exist?

        idx = 1:VAL.feval;
        for g = 1:VAL.n
            idx = idx( MSS.LL(g, idx) - midPoint(g) == 0 );
            if ( isempty( idx ) )
                break
            end
        end

        % Calculate function value in new point and add this point to vertex set
        if ( isempty( idx ) )
            FF = VAL.TPD( midPoint );
            if ( FF <= VAL.Fmin )
                VAL.Fmin = FF;
                VAL.Xmin = midPoint;
            end
            [VAL.feval, idx] = deal( VAL.feval + 1 );
            MSS.LL(:, VAL.feval) = [midPoint; FF];
        else
            FF = MSS.LL(VAL.nn, idx);
        end

        % Create subdivided simplices
        MSS.SS(VAL.I).V = MSS.SS(POH(i)).V;
        [MSS.SS(POH(i)).V(:, kk), MSS.SS(VAL.I).V(:, ll)] = deal( idx );
        [MSS.SS(POH(i)).VV(kk), MSS.SS(VAL.I).VV(ll)] = deal( FF );
        MSS.CC(:, POH(i)) = VAL.alpha( MSS.LL(1:VAL.n, MSS.SS(POH(i)).V) );
        MSS.CC(:, VAL.I) = VAL.alpha( MSS.LL(1:VAL.n, MSS.SS(VAL.I).V) );

        MSS.DD(POH(i)) = VAL.DiamL( POH(i), MSS );
        MSS.DD(VAL.I) = VAL.DiamL( VAL.I, MSS );
        MSS.FF(POH(i)) = VAL.Fvial( POH(i), MSS );
        MSS.FF(VAL.I) = VAL.Fvial( VAL.I, MSS );

        if ( t ~= Dim )

            if ( min( MSS.SS(VAL.I).VV ) < min( MSS.SS(POH(i)).VV ) )
                POH(i) = VAL.I;
            end

            V = MSS.LL(1:VAL.n, MSS.SS(POH(i)).V);
            longDistD1 = roundn( pdist2( V', V' ), -16 );
            [LG1, kk] = max( max( longDistD1, [], 1 ) );

            if ( LG1 == LG )
                [~, ll] = max( longDistD1(:, kk) );
            else
                break;
            end

        end

        % Break if function evaluations exceed the budget
        if ( VAL.Limit <= VAL.I )
            break;
        end

    end

    % Break if function evaluations exceed the budget
    if ( VAL.Limit <= VAL.I )
        break;
    end

end

return