function [ MSS, VAL ] = DTCS( VAL, MSS, POH )

for i = 1:size( POH, 2 )

    Dim = VAL.condition * VAL.n + ~VAL.condition;
    longDistD = roundn( pdist2( MSS.Simpl(POH(i)).V', MSS.Simpl(POH(i)).V' ), -16 );
    [LG, kk] = max( max( longDistD, [], 1 ) );
    [~, ll] = max( longDistD(:, kk) );

    for t = 1:Dim

        % Indexes
        C = POH(i);
        A = VAL.I + 1;
        [B, VAL.I] = deal( VAL.I + 2 );

        % New left vertex point
        m1 = 1/3*MSS.Simpl(C).V(:, kk) + 2/3*MSS.Simpl(C).V(:, ll);

        % New right vertex point
        m2 = 2/3*MSS.Simpl(C).V(:, kk) + 1/3*MSS.Simpl(C).V(:, ll);

        [MSS.Simpl(A).V, MSS.Simpl(B).V] = deal( MSS.Simpl(C).V );
        [MSS.Simpl(C).V(:, ll), MSS.Simpl(B).V(:, kk)] = deal( m1 );
        [MSS.Simpl(C).V(:, kk), MSS.Simpl(A).V(:, ll)] = deal( m2 );

        % Calculate new points
        MSS.CC(:, A) = sum( VAL.alpha * MSS.Simpl(A).V' )';
        MSS.CC(:, B) = sum( VAL.alpha * MSS.Simpl(B).V' )';

        % Transform and evaluate points at the objective function
        MSS.FF(A) = VAL.TPD( MSS.CC(:, A) );
        MSS.FF(B) = VAL.TPD( MSS.CC(:, B) );

        % Calculate diameters
        MSS.DD(A) = VAL.DiamL( A, MSS );
        MSS.DD(B) = VAL.DiamL( B, MSS );
        MSS.DD(C) = VAL.DiamL( C, MSS );

        if ( t ~= Dim )
            
            ii = [A, B, C];
            [~, ind] = min( MSS.FF(ii) );
            POH(i) = ii( ind );

            longDistD1 = roundn( pdist2( MSS.Simpl(POH(i)).V', MSS.Simpl(POH(i)).V' ), -16 );
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

% Update globl values
New_samples = flip( [VAL.Fminindex, VAL.feval + 1:VAL.I] );
[ Fmin, min_index ] = min( MSS.FF(New_samples) );
if VAL.Fmin >= Fmin
    VAL.Fmin = Fmin;
    VAL.Fminindex = New_samples(min_index);
    VAL.Xmin = MSS.CC(:, VAL.Fminindex);
end
VAL.feval = VAL.I;

return