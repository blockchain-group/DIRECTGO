function [ II, Dat ] = MultiLevelSearch( Dat )

% Find good indexes
II = find(Dat.MSS.DD(1:Dat.VAL.I) ~= -1);

% Update Fmin_old
Dat.VAL.Fmin_old = Dat.VAL.Fmin;

return