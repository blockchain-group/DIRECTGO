function [ MSS, VAL ] = DTDV( VAL, MSS, POH )

for j = 1:size( POH, 2 )

    % Find long side
    max_L = max( MSS.LL(:, MSS.RECT(POH(j)).h) );

    % Split sides
    lls = VAL.condition * length( find( MSS.LL(:, MSS.RECT(POH(j)).h) == max_L ) ) + ~VAL.condition;
    for i = 1:lls

        % Indexes
        A = VAL.I + 1;
        [B, VAL.I] = deal( VAL.I + 2 );
        C = POH(j);
        ii = [C, A, B];

        % Find dimension(s) for trisection
        ls = find( MSS.LL(:, MSS.RECT(C).h) == max_L, 1, 'first' );
        set = [ls, 1:(ls - 1), (ls + 1):VAL.n];
        MSS.RECT(C).h = MSS.RECT(C).h + 1;
        delta = max_L / 3;
        MSS.LL(:, MSS.RECT(C).h) = MSS.LL(:, MSS.RECT(C).h - 1);
        MSS.LL(ls, MSS.RECT(C).h) = delta;

        % New sampling points
        CSS1 = repmat( MSS.CX(:, MSS.RECT(C).p1), 1, 3 );
        CSS2 = repmat( MSS.CX(:, MSS.RECT(C).p2), 1, 3 );

        sign = 2 * (CSS1(ls, 1) < CSS2(ls, 1)) - 1;
        CSS1(ls, 2:3) = CSS1(ls, 2:3) + sign * 2 * delta;
        CSS2(ls, 1:2) = CSS2(ls, 1:2) - sign * 2 * delta;

        MSS.RECT([A, B]) = repmat( MSS.RECT(C), 1, 2 );

        left = 1:VAL.feval;
        right = left;
        for g = set
            left = left( MSS.CX(g, left) - CSS1(g, 2) == 0 );
            if ( isempty( left ) )
                break
            end
        end
        [ MSS.RECT(A).p1, MSS.RECT(B).p1 ] = deal( left );

        if ( isempty( MSS.RECT(A).p1 ) )
            [ VAL.feval, MSS.RECT(A).p1, MSS.RECT(B).p1 ] = deal( VAL.feval + 1 );
            MSS.CX(:, VAL.feval) = CSS1(:, 2);
            MSS.FX(VAL.feval) = VAL.TPD( CSS1(:, 2) );
            if ( MSS.FX(VAL.feval) < VAL.Fmin )
                VAL.Fmin = MSS.FX(VAL.feval);
                VAL.Xmin = CSS1(:, 2);
            end
        end

        for g = set
            right = right( MSS.CX(g, right) - CSS2(g, 2) == 0 );
            if ( isempty( right ) )
                break
            end
        end
        [ MSS.RECT(A).p2, MSS.RECT(C).p2 ] = deal( right );

        if ( isempty( MSS.RECT(A).p2 ) )
            [ VAL.feval, MSS.RECT(A).p2, MSS.RECT(C).p2 ] = deal( VAL.feval + 1 );
            MSS.CX(:, VAL.feval) = CSS2(:, 2);
            MSS.FX(VAL.feval) = VAL.TPD( CSS2(:, 2) );
            if ( MSS.FX(VAL.feval) < VAL.Fmin )
                VAL.Fmin = MSS.FX(VAL.feval);
                VAL.Xmin = CSS2(:, 2);
            end
        end

        % Store \ update information
        MSS.FF(C) = VAL.Fvial( MSS.FX(MSS.RECT(C).p1), MSS.FX(MSS.RECT(C).p2) );
        MSS.FF(A) = VAL.Fvial( MSS.FX(MSS.RECT(A).p1), MSS.FX(MSS.RECT(A).p2) );
        MSS.FF(B) = VAL.Fvial( MSS.FX(MSS.RECT(B).p1), MSS.FX(MSS.RECT(B).p2) );
        MSS.DD(ii) = repmat( VAL.DiamL(MSS.LL(:, MSS.RECT(C).h)), 1, 3 );
        MSS.CC(:, ii) = ( CSS1(:, 1:3) + CSS2(:, 1:3) ) / 2;

        [~, ind] = min( MSS.FF(ii) );
        POH(j) = ii( ind );

        % Break if function evaluations exceed the budget
        if ( VAL.Limit <= VAL.feval )
            break;
        end

    end

    % Break if function evaluations exceed the budget
    if ( VAL.Limit <= VAL.feval )
        break;
    end

end

return