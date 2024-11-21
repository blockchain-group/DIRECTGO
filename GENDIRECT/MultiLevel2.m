function [ II, Dat ] = MultiLevel2( Dat )

% Find good indexes
II = find( Dat.MSS.DD(1:Dat.VAL.I) ~= -1 );

% Find W level
W_c = Dat.VAL.Wcicle( Dat.VAL.countG );
W_Level = (W_c == 2) * 1 + (W_c == 1) * (9/10) + (W_c == 0) * (1/9);

% Find indexes at specific W level
II = II( Dat.MSS.DD(II) <= quantile( Dat.MSS.DD(II), W_Level ) );

% Fix ep parameter
Dat.Selection.Ep = (W_Level == 2) * 1e-5 + (W_Level == 1) * 1e-7;

% Move to next level
Dat.VAL.countG = mod( Dat.VAL.countG, 8 ) + 1;

% Update Fmin_old
Dat.VAL.Fmin_old = Dat.VAL.Fmin;

return