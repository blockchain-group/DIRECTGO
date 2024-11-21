function Dat = Constructor(Dat)

    switch Dat.Partitioning.Strategy
        case 'DBC'
            Dat.SubDivide = @DBC;
            Dat.Initialization = @DBC_Init;
            Dat.STA.Dbc = true;
        case 'DBVD'
            Dat.SubDivide = @DBVD;
            Dat.Initialization = @DBVD_Init;
        case 'DTCS'
            Dat.SubDivide = @DTCS;
            Dat.Initialization = @DTCS_Init;
        case 'DBDP'
            Dat.SubDivide = @DBDP;
            Dat.Initialization = @DBDP_Init;
        case 'DTDV'
            Dat.SubDivide = @DTDV;
            Dat.Initialization = @DTDV_Init;
        case 'DBVS'
            Dat.SubDivide = @DBVS;
            Dat.Initialization = @DBVS_Init;
        case 'DTC'
            Dat.SubDivide = @DTC;
            Dat.Initialization = @DTC_Init;
        otherwise
            fprintf('No such option ''%s'' for Partitioning.Strategy. Please choose one of the following options: ''DBC'', ''DBVS'', ''DTCS'', ''DBDP'', ''DTDV'', ''DBDPV'', ''DTC'', \n', Dat.Partitioning.Strategy);
                Dat.perror = -1;
    end
    switch Dat.Selection.Strategy
        case 'Aggressive'
            Dat.Selection.Strategy = @Find_Aggressive_POH;
            Dat.STA.Ags = true;
        case 'Pareto'
            Dat.Selection.Strategy = @Find_Pareto_POH;
        case 'RedPareto'
            Dat.Selection.Strategy = @Find_RedPareto_POH;
        case 'Original'
            Dat.Selection.Strategy = @Find_Orinial_POH;
        otherwise
            fprintf('No such option ''%s'' for Selection.Strategy. Please choose one of the following options: ''Aggressive'', ''Pareto'', ''RedPareto'', ''Original'', \n', Dat.Selection.Strategy);
            Dat.perror = -1;
    end
    if ~strcmp('Off', Dat.Selection.ControlEp)
        switch Dat.Selection.ControlEp
            case 'Restart'
                Dat.Selection.ControlEp = @Restart;
            case 'MultiLevel1'
                Dat.Selection.ControlEp = @MultiLevel1;
            case 'MultiLevel2'
                Dat.Selection.ControlEp = @MultiLevel2;
            otherwise
                fprintf('No such option ''%s'' for Dat.Selection.ControlEp. Please choose one of the following options: ''Restart'', ''MultiLevel1'', ''MultiLeve2'', \n', Dat.Selection.ControlEp);
                Dat.perror = -1;
        end
    else
        Dat.Selection.ControlEp = @MultiLevelSearch;
    end
    if ~strcmp('Off', Dat.Selection.SolRefin)
        switch Dat.Selection.SolRefin
            case 'Min'
                Dat.Selection.SolRefin = @LimitMin;
            case 'Median'
                Dat.Selection.SolRefin = @LimitMedian;
            case 'Average'
                Dat.Selection.SolRefin = @LimitAverage;
            otherwise
                fprintf('No such option ''%s'' for Dat.Selection.SolRefin. Please choose one of the following options: ''Min'', ''Median'', ''Average'', \n', Dat.Selection.SolRefin);
                Dat.perror = -1;
        end
    else
        Dat.Selection.SolRefin = @Refinement;
    end
    if ~strcmp('Off', Dat.Hybridization.Strategy)
        Dat.Hybridization.options = optimoptions('fmincon', 'Display', 'none', 'MaxFunctionEvaluations', Dat.Hybridization.MaxEvaluations, 'MaxIterations', Dat.Hybridization.MaxIterations, 'Algorithm', Dat.Hybridization.LocalSearch);
        switch Dat.Hybridization.Strategy
            case 'Aggressive'
                Dat.Hybridization.x0 = [];
                Dat.Hybridization.Strategy = @Aggressive;
            case 'Single'
                Dat.Hybridization.Strategy = @Single;
            case 'Clustering'
                Dat.Hybridization.Strategy = @Clustering;
            otherwise
                fprintf('No such option ''%s'' for Hybridization.Strategy. Please choose one of the following options: ''Aggressive'', ''Single'', ''Clustering'', ''Off'', \n', Dat.Hybridization.Strategy);
                Dat.perror = -1;
        end
    else
        Dat.Hybridization.Strategy = @(x) x;
    end
    if strcmp('On', Dat.Selection.GloballyBiased)
        Dat.STA.Gb = true;
    end
    if strcmp('All', Dat.Selection.EqualCand)
        Dat.STA.Eqs = true;
    end
    if strcmp('On', Dat.Selection.TwoPhase)
        Dat.STA.TPgl = true;
    end
end