function [ MSS, VAL ] = DBC( VAL, MSS, POH )

for g = 1:size( POH, 2 )

    % Find long side
    max_L = max( MSS.LL(:, POH(g)) );

    % Split sides
    lls = VAL.condition * length( find( MSS.LL(:, POH(g)) == max_L ) ) + ~VAL.condition;
    for i = 1:lls

        % Index of the center point and new ones
        L = VAL.I + 1;
        R = POH(g);
        POHd = MSS.LFA(1, R);
        id = VAL.feval + 1:VAL.feval + 2;

        % Find dimension(s) for trisection
        ls = find( MSS.LL(:, R) == max_L );

        if ( size( ls, 1 ) ~= 1 )

            % Find closest to the Xmin
            DIR = max( abs( MSS.CC(:, VAL.Fminindex) - MSS.CC(:, POHd) ), [], 2 );
            ls  = ls(find( DIR(ls) == max(DIR(ls) ), 1, 'first'));

        end

        % Calculate side lengths
        MSS.LL(ls, R) = max_L / 2;
        MSS.LL(:, L) = MSS.LL(:, R);

        % Calculate new points
        MSS.CC(:, id) = repmat(MSS.CC(:, POHd), 1, 2);
        MSS.CC(ls, id(1)) = MSS.CC(ls, id(1)) - max_L / 4;
        MSS.CC(ls, id(2)) = MSS.CC(ls, id(2)) + max_L / 4;

        % Evaluate the objective function
        MSS.FE(id(1)) = VAL.TPD( MSS.CC(:, id(1)) );
        MSS.FE(id(2)) = VAL.TPD( MSS.CC(:, id(2)) );

        % Update Fmin and Xmin
        ids = flip( [VAL.Fminindex, id] );
        [ Fmin, min_index ] = min( MSS.FE(ids) );
        if VAL.Fmin >= Fmin
            VAL.Fmin = Fmin;
            VAL.Fminindex = ids(min_index);
            VAL.Xmin = MSS.CC(:, VAL.Fminindex);
        end

        % Calculate diameters
        MSS.DD([ L, R ]) = VAL.DiamL( MSS.LL(:, L) );

        % Update approximation points
        HA = MSS.LFA(2:(find( isnan( MSS.LFA(:, R) ), 1, 'first' ) - 1), R);
        l_HA = [ id(1); HA( MSS.CC(ls, HA) <= MSS.CC(ls, POHd) ); POHd ];
        r_HA = [ id(2); HA( MSS.CC(ls, HA) >= MSS.CC(ls, POHd) ); POHd ];

        % Store \ update information
        MSS.LFA(1:length(l_HA), R) = l_HA;
        MSS.LFA(1:length(r_HA), L) = r_HA;
        MSS.LFmin(1:length(l_HA), R) = sort( MSS.FE(l_HA) )';
        MSS.LFmin(1:length(r_HA), L) = sort( MSS.FE(r_HA) )';
        MSS.FF(R) = VAL.Fvial( MSS.FE(l_HA), MSS.FE(id(1)) );
        MSS.FF(L) = VAL.Fvial( MSS.FE(r_HA), MSS.FE(id(2)) );

        set = [ L, R ];
        [ ~, ind ] = min( MSS.FE(set) );
        POH(g) = set(ind);
        VAL.I = L;
        VAL.feval = VAL.feval + 2;

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