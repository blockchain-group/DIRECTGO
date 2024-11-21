function [ MSS, VAL ] = DBDP( VAL, MSS, POH )

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
        MSS.RECT(R).h = MSS.RECT(R).h + 1;
        VAL.e = max( [ VAL.e, MSS.RECT(R).h ] );
        if ( VAL.e == MSS.RECT(R).h )
            MSS.LL(:, MSS.RECT(R).h)  = MSS.LL(:, MSS.RECT(R).h - 1);
            MSS.LL(ls, MSS.RECT(R).h) = max_L / 2;
        end
        MSS.RECT(L) = MSS.RECT(R);

        % New sampling points
        sign = 2 * (MSS.RECT(R).p1(ls) < MSS.RECT(R).p2(ls)) - 1;
        MSS.RECT(L).p1(ls) = MSS.RECT(R).p1(ls) + sign * max_L / 2;
        MSS.RECT(R).p2(ls) = MSS.RECT(L).p2(ls) - sign * max_L / 2;

        % Transform back to the original space and evaluate the function
        MSS.RECT(L).f1 = VAL.TPD( MSS.RECT(L).p1 );
        MSS.RECT(R).f2 = VAL.TPD( MSS.RECT(R).p2 );
        MSS.FF(R) = VAL.Fvial( MSS.RECT(R).f1, MSS.RECT(R).f2 );
        MSS.FF(L) = VAL.Fvial( MSS.RECT(L).f1, MSS.RECT(L).f2 );

        % Update Fmin and Xmin
        set = [ MSS.RECT(R).f1, MSS.RECT(R).f2, MSS.RECT(L).f1, MSS.RECT(L).f2 ];
        [ vv, ii ] = min( set );
        if ( min( set ) < VAL.Fmin )
            setc = [ MSS.RECT(R).p1, MSS.RECT(R).p2, MSS.RECT(L).p1, MSS.RECT(L).p2 ];
            VAL.Fmin = vv;
            VAL.Xmin = setc(:, ii);
        end

        % Store \ update information
        MSS.CC(:, L) = ( MSS.RECT(L).p1 + MSS.RECT(L).p2 ) / 2;
        MSS.CC(:, R) = ( MSS.RECT(R).p1 + MSS.RECT(R).p2 ) / 2;
        MSS.DD([ R, L ]) = deal( VAL.DiamL( MSS.LL(:, MSS.RECT(R).h) ) );

        [ ~, ind ] = min( [ MSS.FF(R), MSS.FF(L) ] );
        io = [ R, L ];
        POH(j) = io(ind);
        VAL.feval = VAL.feval + 2;
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