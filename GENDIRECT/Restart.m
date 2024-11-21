function [ II, Dat ] = Restart( Dat )

% Find good indexes
II = find( Dat.MSS.DD(1:Dat.VAL.I) ~= -1 );

% Solution stagnated?
CC = ( abs( Dat.VAL.Fmin_old - Dat.VAL.Fmin ) < 10^(-4) );
Dat.VAL.Fmin_stag = Dat.VAL.Fmin_stag * CC + CC;

% Fix ep parameter
Dat.Selection.Ep = (Dat.VAL.Fmin_stag > 5 && Dat.VAL.Fmin_stag < 50) * 1e-2;

% Update Fmin_old
Dat.VAL.Fmin_old = Dat.VAL.Fmin;

return