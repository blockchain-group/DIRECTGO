function [ MSS, VAL ] = DBVD( VAL, MSS, POH )

for j = 1:size( POH, 2 )
    % Find long side
    max_L = max( MSS.LL(:, MSS.RECT(POH(j)).h) );

    % Split sides
    lls = VAL.condition * length( find( MSS.LL(:, MSS.RECT(POH(j)).h) == max_L ) ) + ~VAL.condition;
    for i = 1:lls

        % Left and right indexes
        L = VAL.I + 1;
        R = POH(j);

        % Find dimension(s) for trisection
        ls = find( MSS.LL(:, MSS.RECT(R).h) == max_L, 1, 'first' );
        set = [ls, 1:(ls - 1), (ls + 1):VAL.n];
        delta = max_L / 2;
        MSS.RECT(R).h = MSS.RECT(R).h + 1;
        MSS.LL(:, MSS.RECT(R).h) = MSS.LL(:, MSS.RECT(R).h - 1);
        MSS.LL(ls, MSS.RECT(R).h) = delta;
        MSS.RECT(L) = MSS.RECT(R);

        % New sampling points
        MSS.CC(:, L) = MSS.CC(:, R);
        [ vertex, vertex2 ] = deal( MSS.CX(:, MSS.RECT(R).p2) );
        sign = 2 * (MSS.RECT(R).p1(ls) < vertex(ls)) - 1;
        MSS.CC(ls, R) = MSS.CC(ls, R) + sign * delta / 2;
        MSS.CC(ls, L) = MSS.CC(ls, L) - sign * delta / 2;
        MSS.RECT(L).p1(ls) = MSS.RECT(R).p1(ls) + sign * delta * (2 / 3);
        vertex(ls) = vertex(ls) - sign * delta * 2;

        MSS.RECT(R).p2 = 1:VAL.e;
        for g = set
            MSS.RECT(R).p2 = MSS.RECT(R).p2( MSS.CX(g, MSS.RECT(R).p2) - vertex(g) == 0 );
            if ( isempty( MSS.RECT(R).p2 ) )
                break
            end
        end

        if ( isempty( MSS.RECT(R).p2 ) )
            VAL.e = VAL.e + 1;
            VAL.feval = VAL.feval + 1;
            MSS.FX(VAL.e) = VAL.TPD( vertex );
            MSS.CX(:, VAL.e) = vertex;
            MSS.RECT(R).p2 = VAL.e;
        end

        MSS.RECT(L).f1 = VAL.TPD( MSS.RECT(L).p1 );
        MSS.FF(R) = VAL.Fvial( MSS.RECT(R).f1, MSS.FX(MSS.RECT(R).p2) );
        MSS.FF(L) = VAL.Fvial( MSS.RECT(L).f1, MSS.FX(MSS.RECT(L).p2) );

        % Update Fmin and Xmin
        set = [ MSS.RECT(R).f1, MSS.FX(MSS.RECT(R).p2), MSS.RECT(L).f1, MSS.FX(MSS.RECT(L).p2) ];
        [ vv, ii ] = min( set );
        if ( vv < VAL.Fmin )
            setc = [ MSS.RECT(R).p1, vertex, MSS.RECT(L).p1, vertex2 ];
            VAL.Fmin = vv;
            VAL.Xmin = setc(:, ii);
        end

        % Store \ update information
        MSS.DD([ R, L ]) = deal( VAL.DiamL( MSS.LL(:, MSS.RECT(R).h) ) );
        [ ~, ind ] = min( [ MSS.FF(R), MSS.FF(L) ] );
        io = [R, L];
        POH(j) = io(ind);
        VAL.feval = VAL.feval + 1;
        VAL.I = VAL.I + 1;

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