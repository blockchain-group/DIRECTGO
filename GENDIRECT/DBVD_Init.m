function [ MSS, VAL ] = DBVD_Init( VAL, Selection, Partitioning )

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
        VAL.DiamL = @(x) round( (2/3) * norm( x, 2 ), 16 );
end

% Aggregated function value
switch Partitioning.AggrFuncVal
    case 'Mean'
        VAL.Fvial = @(x, y) mean( [ x, y ] );
    otherwise
        VAL.Fvial = @(x, y) min( [ x, y ] );
end

% Initialize
[VAL.I, VAL.e]  = deal( 1 );
VAL.feval = 2;
MSS.RECT(1).p1 = ones( VAL.n, 1 ) ./ 3;
MSS.RECT(1).p2 = 1;
MSS.RECT(1).h  = 1;
MSS.RECT(1).f1 = VAL.TPD( MSS.RECT(1).p1 );
MSS.FX(1) = VAL.TPD( ones( VAL.n, 1 ) );

if ( MSS.RECT(1).f1 < MSS.FX(1) )
    [ VAL.Fmin, VAL.Fmin_old ] = deal( MSS.RECT(1).f1 );
    VAL.Xmin = MSS.RECT(1).p1;
else
    [ VAL.Fmin, VAL.Fmin_old ] = deal( MSS.FX(1) );
    VAL.Xmin = MSS.RECT(1).p2;
end

MSS.FF(1) = VAL.Fvial( MSS.RECT(1).f1, MSS.FX(1) );
MSS.CX(:, 1) = ones(VAL.n, 1);
MSS.CC(:, 1) = ones(VAL.n, 1)./2;
MSS.DD(1) = VAL.DiamL( MSS.LL(:, 1) );

return