function [ MSS, VAL ] = DBVS_Init( VAL, Selection, Partitioning )

% Allocate space
VAL.nn = VAL.n + 1;
VAL.alpha = @(x) sum( (1 / VAL.nn) * x' )';

MSS = struct( 'FF',  nan( 1, VAL.Limit ), ...
    'DD', -ones( 1, VAL.Limit ), ...
    'LL',  zeros( VAL.nn, VAL.Limit ), ...
    'CC',  zeros( VAL.n, VAL.Limit ) ...
    );

% Initialize
V = VertexGeneration( VAL.n );
Simpl = VertexTriangulation( VAL.n )';
VAL.I = factorial( VAL.n );
if ( VAL.I == 1 )
    VAL.feval = VAL.nn;
else
    VAL.feval = length( V );
end
VV = zeros( 1, length( V ) );

% Calculate function values on all vertices V
for i = 1:VAL.feval

    % Transform back to the original space
    VV(i) = VAL.TPD( V(:, i) ); 

    if ( VAL.Fmin > VV(i) )

        [VAL.Fmin_old, VAL.Fmin] = deal( VV(i) );
        VAL.Xmin = V(:, i);

    end

end

for j = 1:VAL.I
    for k = 1:VAL.nn
        for i = 1:VAL.feval
            if ( norm( Simpl(j).V(:, k) - V(:, i), 2 ) == 0 )
                MSS.SS(j).V(k) = i;
                MSS.SS(j).VV(k) = VV(i);
                break;
            end
        end
    end
end

% Find vertex with minimal function value
MSS.LL(:, 1:length(V)) = [V; VV];

for i = 1:VAL.I
    MSS.CC(:, i) = VAL.alpha( Simpl(i).V );
    switch Selection.CandMeasure
        case 'LongSide'
            VAL.DiamL = @(i, MSS) roundn( max( max( pdist2( MSS.LL(1:VAL.n, MSS.SS(i).V)', MSS.LL(1:VAL.n, MSS.SS(i).V)' ) ) ), -16 );
            MSS.DD(i) = VAL.DiamL( i, MSS );
        otherwise
            VAL.DiamL = @(i, MSS) roundn( max( sqrt( sum( (MSS.CC(:, i) - MSS.LL(1:VAL.n, MSS.SS(i).V)) .^ 2, 1 ) ) ), -16 );
            MSS.DD(i) = VAL.DiamL( i, MSS );
    end

    switch Partitioning.AggrFuncVal
        case 'Mean'
            VAL.Fvial = @(i, MSS) mean( MSS.SS(i).VV );
            MSS.FF(i) = VAL.Fvial( i, MSS );
        otherwise
            VAL.Fvial = @(i, MSS) min( MSS.SS(i).VV );
            MSS.FF(i) = VAL.Fvial( i, MSS );
    end
end

return