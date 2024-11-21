function Dat = GloballyBasedstr(Dat)

    condition = strcmp('On', Dat.Selection.GloballyBiased);
    condition = condition & (abs(Dat.VAL.Fmin - Dat.VAL.fMinBeforeImpr) > 0.01 * abs(Dat.VAL.Fmin));
    Dat.VAL.fMinNotImpr(condition) = 0;
    Dat.VAL.fMinBeforeImpr(condition) = Dat.VAL.Fmin;

end