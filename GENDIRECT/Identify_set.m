function [ SET, Dat ] = Identify_set( Dat )

% Prepare for the POHs selection
[ II, Dat ] = Dat.Selection.ControlEp( Dat );

% Find first POH on the plane
[ i_min, Dat, C ] = Dat.Selection.SolRefin(Dat.MSS.FF(II), Dat.MSS.DD(II), true, Dat);

% Find the set of POHs
SET = II(Dat.Selection.Strategy(Dat.MSS.FF(II), C, II, Dat, i_min));

if ( Dat.STA.TPgl )

    % Find Local set of POHs
    Main_set = sum( (Dat.VAL.Xmin - Dat.MSS.CC(:, II) ).^2, 1).^0.5;
    [ i_min, Dat, C ] = Dat.Selection.SolRefin(Main_set, Dat.MSS.DD(II), false, Dat);
    SET2 = II(Dat.Selection.Strategy(Main_set, C, II, Dat, i_min));
    SET = unique( [ SET, SET2 ] );
    SET = sortrows( [ SET; Dat.MSS.DD(SET) ].', 2 ).';

end

end

