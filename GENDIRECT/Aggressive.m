function Dat = Aggressive(Dat)

% Run local search from each POH
[xminlocs,fminlocs,~,outputs]=arrayfun(@(idx)fmincon(@(x)Dat.VAL.TPD(x),...
    Dat.MSS.CC(:,Dat.POHLocal(idx)),[],[],[],[],zeros(Dat.Problem.n,1), ...
    ones(Dat.Problem.n,1),[],Dat.Hybridization.options), ...
    1:length(Dat.POHLocal),'UniformOutput',false);

% Extract data from local solver
Dat.VAL.nLocSearch = Dat.VAL.nLocSearch + length(Dat.POHLocal);
Dat.VAL.fevalLocal = Dat.VAL.fevalLocal + sum(cellfun(@(out) out.funcCount, outputs));
[mins, idxs] = min(cell2mat(fminlocs));

%  Update Fmin if it has improved
if mins < Dat.VAL.Fmin
    Dat.VAL.Fmin = mins;
    Dat.VAL.Xmin = (xminlocs{idxs} - Dat.Problem.xl)./Dat.VAL.delta;
end

end