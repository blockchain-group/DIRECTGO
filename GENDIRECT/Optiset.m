function [VAL, MSS, Problem, optParam] = Optiset(Problem, VAL, optParam, Division)

    tic
    if Problem.info == true
        getInfo = feval(Problem.f);
        if isfield(getInfo, 'confun'), Problem.constraint = getInfo.confun; end
        if getInfo.nx ~= 0, Problem.n = getInfo.nx; end
        Problem.xl = getInfo.xl(Problem.n);
        Problem.xu = getInfo.xu(Problem.n);
        globalMIN  = getInfo.fmin(Problem.n);
        if globalMIN ~= 0
            optParam.goal = globalMIN + abs(globalMIN*optParam.tol)/100;
        else
            optParam.goal = optParam.tol/100;
        end
    else
        Problem.xl = Problem.xL;
        Problem.xu = Problem.xU;
        optParam.goal = Problem.fgoal;
    end
    VAL.PF = 'Iteration: %4i  min value: %15.10f    time(s): %10.05f  Local runs: %3i function evaluations: %8i\n';
    VAL.delta = abs(Problem.xu - Problem.xl);
    VAL.Con = inf;
    VAL.perror = 10;
    VAL.Limit = optParam.maxevals;
    VAL.limits = 1000;
    VAL.condition = strcmp(Division.SubSides, 'All');
    VAL.POHLocal = [];
    VAL.GLC = Problem.n*100 + 1;
    VAL.GLCcount = 0;
    VAL.fevalLocal = 0;
    VAL.nLocSearch  = 0;
    VAL.ep = 0;
    VAL.z =(2*Problem.n) + 1; 
    VAL.fMinBeforeImpr = inf;
    VAL.fMinNotImpr = 0;
    VAL.Fmin_stag = 0;
    VAL.time = 0;
    VAL.z =(2*Problem.n) + 1;
    VAL.Fmin_old = inf;
    VAL.Fmin    = inf;
    VAL.I = 1;
    VAL.n = Problem.n;
    VAL.Xmin    = zeros(Problem.n, 1);
    VAL.minas   = zeros(Problem.n, 1);
    VAL.TPD = @(x) Evaluate( Problem, bsxfun( @times, VAL.delta, x ) + Problem.xl );
    MSS = struct('FF', nan, 'XX', nan, 'df', nan, 'LFmin', nan(2*Problem.n + 1, 1), 'DD', -1, 'CC', zeros(Problem.n, 1)); 

return