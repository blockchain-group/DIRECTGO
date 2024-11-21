function Dat = Single(Dat)

% Check if DIRECT sufficiently improved Fmin
if (abs(Dat.VAL.Fmin - Dat.VAL.Fmin_old) > 0.01*abs(Dat.VAL.Fmin)) || Dat.VAL.itctr == -1

    % Run local search from best point
    [xminloc,fminloc,~,output]=fmincon(@(x)Dat.VAL.TPD(x),Dat.VAL.Xmin,...
        [],[],[],[],zeros(Dat.Problem.n,1),ones(Dat.Problem.n,1),[],...
        Dat.Hybridization.options);

    %  Update Fmin if it has improved
    Dat.VAL.nLocSearch = Dat.VAL.nLocSearch + 1;
    Dat.VAL.fevalLocal = Dat.VAL.fevalLocal + output.funcCount;
    if fminloc < Dat.VAL.Fmin
        Dat.VAL.Fmin = fminloc;
        Dat.VAL.Xmin = (xminloc - Dat.Problem.xl)./Dat.VAL.delta;
    end

end

end

