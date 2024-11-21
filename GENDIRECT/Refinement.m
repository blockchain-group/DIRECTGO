function [i_min, Dat, d] = Refinement(fc, szes, Flag, Dat)

% Find the first POH 
[ ~, i_min ] = min(fc);
d = unique( szes );

if ( Dat.STA.Ags )
    idx = 1;
else
    idx = find( d == szes( i_min ) );
end

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