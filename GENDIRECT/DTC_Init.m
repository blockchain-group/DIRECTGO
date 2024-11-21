function [ MSS, VAL ] = DTC_Init( VAL, Selection, ~ )

% Allocate space
MSS = struct( 'FF',  nan( 1, VAL.Limit ), ...
    'DD', -ones( 1, VAL.Limit ), ...
    'LL',  ones( VAL.n, VAL.Limit ), ...
    'CC',  zeros( VAL.n, VAL.Limit ) ...
    );

% Hyper-rectangle size measure
switch Selection.CandMeasure
    case 'LongSide'
        VAL.DiamL = @(x) round( norm( x, "inf" ), 16 );
    otherwise
        VAL.DiamL = @(x) round( 0.5 * norm( x, 2 ), 16 );
end

% Initialize
[ VAL.I, VAL.feval, VAL.Fminindex ] = deal( 1 );
MSS.DD(1) = VAL.DiamL( MSS.LL(:, 1) );
[ MSS.CC(:, 1), VAL.Xmin ] = deal( ones(VAL.n, 1) / 2 );
MSS.FF(1) = VAL.TPD( MSS.CC(:, 1) );
[ VAL.Fmin_old, VAL.Fmin ] = deal( MSS.FF(1) );

return