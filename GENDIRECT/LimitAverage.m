function [i_min, Dat, d] = LimitAverage(fc, szes, Flag, Dat)

% Find the first POH
if ( Flag )
    E = max( Dat.Selection.Ep * abs( Dat.VAL.Fmin - mean( fc ) ), 1e-8 );
else
    E = max( Dat.Selection.Ep * abs( mean( fc ) ), 1e-8 );
end
[ ~, i_min ] = min( (fc - min( fc ) + E )./ szes );
d = unique( szes );
idx = find( d == szes(i_min), 1 );

% GloballyBiased strategy
if ( ( Flag ) && ( Dat.STA.Gb ) )
    if ( abs( Dat.VAL.Fmin - fc(i_min) ) < 1e-2 )
        condition = ( abs( Dat.VAL.Fmin - Dat.VAL.fMinBeforeImpr ) <= 1e-2 * abs( Dat.VAL.Fmin ) );
        Dat.VAL.fMinNotImpr = Dat.VAL.fMinNotImpr + condition;
        Dat.VAL.fMinNotImpr = Dat.VAL.fMinNotImpr * condition;
        Dat.VAL.fMinBeforeImpr( condition ) = Dat.VAL.Fmin;
    end
end

if ( Dat.STA.Gb )
    if ( ( Dat.VAL.fMinNotImpr > 5 ) && ( mod( Dat.VAL.fMinNotImpr, 5 ) ~= 0 ) )
        idx = idx + round( 6 * (length( d ) - idx) / 8 );
    end
end
d = d( idx:end );

return