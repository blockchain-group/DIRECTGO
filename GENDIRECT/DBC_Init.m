function [ MSS, VAL ] = DBC_Init( VAL, Selection, Partitioning )

% Allocate space
Fixed_Limit = fix( VAL.Limit / 2 );
MSS = struct( 'FE',    nan( 1, VAL.Limit ), ...
    'FF',    nan( 1, Fixed_Limit ), ...
    'DD',   -ones( 1, Fixed_Limit ), ...
    'LL',    ones( VAL.n, Fixed_Limit ), ...
    'CC',    zeros( VAL.n, VAL.Limit ), ...
    'LFA',   nan( (2*VAL.n) + 1, Fixed_Limit ), ...
    'LFmin', nan( (2*VAL.n) + 1, Fixed_Limit ) ...
    );

% Hyper-rectangle size measure
switch Selection.CandMeasure
    case 'LongSide'
        VAL.DiamL = @(x) round( norm( x, "inf" ), 16 );
    otherwise
        VAL.DiamL = @(x) round( norm( x, 2 ), 16 );
end

% Aggregated function value
switch Partitioning.AggrFuncVal
    case 'Midpoint'
        VAL.Fvial = @(x, y) y;
    case 'Mean'
        VAL.Fvial = @(x, y) mean( x );
    case 'Minimum'
        VAL.Fvial = @(x, y) min( x );
    otherwise
        VAL.Fvial = @(x, y) mean( [ min( x ), y ] );
end

% Initialize
[ VAL.I, VAL.feval, VAL.Fminindex, MSS.LFA(1, 1) ] = deal( 1 );
[ MSS.CC(:, 1), VAL.Xmin ] = deal( ones( VAL.n, 1 ) / 2 );
MSS.DD(1) = VAL.DiamL( MSS.LL(:, 1) );
[ MSS.FE(1), MSS.FF(1), VAL.Fmin_old, VAL.Fmin ] = deal( VAL.TPD( MSS.CC(:, 1) ) );

return