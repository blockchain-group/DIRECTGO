function [ MSS, VAL ] = DTDV_Init( VAL, Selection, Partitioning )

% Allocate space
MSS = struct( 'FF',  nan( 1, VAL.Limit ), ...
    'DD', -ones( 1, VAL.Limit ), ...
    'LL',  ones( VAL.n, VAL.Limit ), ...
    'CC',  zeros( VAL.n, VAL.Limit ), ...
    'CX', -ones( VAL.n, VAL.Limit ), ...
    'FX',  zeros( 1, VAL.Limit ) ...
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
    case 'Mean'
        VAL.Fvial = @(x, y) min( [ x, y ] );
    otherwise
        VAL.Fvial = @(x, y) mean( [ x, y ] );
end

% Initialize
[VAL.I, VAL.e]  = deal( 1 );
VAL.feval = 2;

MSS.RECT(1).h = 1;
MSS.RECT(1).p1 = 1;
MSS.RECT(1).p2 = 2;
MSS.FX(1) = VAL.TPD( zeros( VAL.n, 1 ) );
MSS.FX(2) = VAL.TPD( ones( VAL.n, 1 ) );

MSS.CX(:, 1) = zeros( VAL.n, 1 );
MSS.CX(:, 2) = ones( VAL.n, 1 );
MSS.CC(:, 1) = ( MSS.CX(:, 1) + MSS.CX(:, 2) ) / 2;
MSS.DD(1) = VAL.DiamL( MSS.LL(:, 1) );
MSS.FF(1) = VAL.Fvial( MSS.FX(1), MSS.FX(2) );

[ min_value, min_index ] = min( MSS.FX(1:2) );
[ VAL.Fmin, VAL.Fmin_old ] = deal( min_value );
VAL.Xmin = MSS.CX(:, min_index);

return