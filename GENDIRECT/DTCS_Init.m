function [ MSS, VAL ] = DTCS_Init( VAL, Selection, ~ )

    % Allocate space
    MSS = struct( 'FF',  nan( 1, VAL.Limit ), ...
                  'DD', -ones( 1, VAL.Limit ), ...
                  'CC',  zeros( VAL.n, VAL.Limit ) ...
                 );

    % Initialize
    VAL.nn = VAL.n + 1;
    VAL.alpha = 1 / (VAL.nn);
    MSS.Simpl = VertexTriangulation( VAL.n )';
    VAL.I = factorial( VAL.n );

    for i = 1:VAL.I

        MSS.CC(:, i) = sum( VAL.alpha * MSS.Simpl(i).V' )';
        MSS.FF(i) = VAL.TPD( MSS.CC(:, i) );
        VAL.feval = VAL.feval + 1;

        switch Selection.CandMeasure
            case 'LongSide'
                VAL.DiamL = @(i, MSS) roundn( max( max( pdist2( MSS.Simpl(i).V(:, 1:VAL.nn)', MSS.Simpl(i).V(:, 1:VAL.nn)' ) ) ), -16 );
                MSS.DD(i) = VAL.DiamL( i, MSS );
            otherwise
                VAL.DiamL = @(i, MSS) roundn( max( sqrt( sum( (MSS.CC(:, i) - MSS.Simpl(i).V(:, 1:VAL.nn)) .^ 2, 1 ) ) ), -16 );
                MSS.DD(i) = VAL.DiamL( i, MSS );
        end

    end
    
    % Find min
    [VAL.Fmin, VAL.Fminindex] = min( MSS.FF(1:VAL.I) );
    VAL.Fmin_old = VAL.Fmin;
    VAL.Xmin = MSS.CC(:, VAL.Fminindex);

return