classdef GENDIRECT
    
    properties
        % Full Name of the Algorithm
        name = 'Generalized DIRECT algorithm';
        
        optParam = struct('maxits', 100, 'maxevals', 10000, 'showits', 1, 'goal', -inf, 'tol', 1);
        Partitioning = struct('Strategy', 'DTC', 'SubSides', 'All', 'AggrFuncVal', 'Midpoint'); 
        Selection = struct('Strategy', 'Original', 'CandMeasure', 'Diagonal', 'SolRefin', 'Off', 'EqualCand', 'All', 'Ep', 1e-4, 'Deep', 1e-16, 'ControlEp', 'None', 'GloballyBiased', 'Off', 'TwoPhase', 'Off');
        LocalSelection = struct('Strategy', 'Original', 'SolRefin', 'false', 'EqualCand', 'All', 'Ep', 1e-4, 'Deep', 1e-16, 'ControlEp', 'None', 'GloballyBiased', 'false');
        Hybridization = struct('Strategy', 'Off', 'MaxEvaluations', 3000, 'MaxIterations', 1000, 'LocalSearch', 'sqp'); 
        VAL = struct('time', 0, 'itctr', -1, 'I', 0, 'feval', 0, 'countL', 1, 'countG', 1, 'Wcicle', [2, 1, 0, 1, 1, 0, 1, 2]);
        STA = struct('Ags', false, 'Gb', false, 'Eqs', false, 'Dbc', false, 'TPgl', false);
        Problem = [];
        opts = [];
        perror = 10;
        Initialization = @DTC_Init;
        SubDivide = @DTC;
        
    end

    properties(SetAccess = protected)
        
        Fmin = 0;
        Xmin = [];
        Iterations = 0;
        Time = 0;
        History = [];
        Evaluations = 0;
        POH = [];
        POHLocal = [];
        MSS = [];

    end
    
    properties(Dependent = true)
        
        % Full Name
        full_name = [];

        % Best Objective Value Ever Found
        best_obj_value = [];
        
        % Effective Problem
        eff_problem = [];
        
    end
    
    methods

        % Solves a Optimization Problem
        function Dat = solve(Dat, ~)
            Dat = Constructor(Dat);

            [Dat.VAL, Dat.MSS, Dat.Problem, Dat.optParam] = Optiset(Dat.Problem, Dat.VAL, Dat.optParam, Dat.Partitioning);

            % Global vals
            global Item_History Item_Evals Item_Mmax Item_xMin Item_Iter %#ok<*GVMIS>
            Item_History = [0, 0, inf, 0];
            Item_Evals = 0;
            Item_Iter = 0;
            Item_xMin = inf(Dat.Problem.n, 1);
            Item_Mmax = Dat.optParam.maxevals;

            % Initialization step
            [ Dat.MSS, Dat.VAL ] = Dat.Initialization( Dat.VAL, Dat.Selection, Dat.Partitioning );

            while Dat.VAL.itctr < Dat.optParam.maxits && Dat.VAL.feval + Dat.VAL.fevalLocal < Dat.optParam.maxevals && Dat.VAL.Fmin > Dat.optParam.goal && Dat.perror > 0 && toc < 3600

                % Selection of potential optimal hyper-rectangles step
                Dat = Dat.POH_Selection();

                % Subdivide potential optimal hyper-rectangles
                [ Dat.MSS, Dat.VAL ] = Dat.SubDivide( Dat.VAL, Dat.MSS, Dat.POH );

                Dat = GloballyBasedstr( Dat );

                % Hybridizationize
                Dat = Dat.Hybridization.Strategy(Dat);
                
                % Display
                Dat = Dat.Display();
            end                                                       

            % Return value
            Dat.Stop();
            Dat.Fmin = Item_History(end, 3);
            Dat.Xmin = Item_xMin;
            Dat.Iterations = Item_Iter;
            Dat.Time = Item_History(end, 4);
            Dat.History = Item_History(2:end, :);
            Dat.Evaluations = Item_Evals;
        end
        
    end

    methods(Access = protected)
        
        function Dat = Display(Dat)
            global Item_Iter Item_History
            [Dat.VAL.itctr, Item_Iter] = deal(Dat.VAL.itctr + 1);
            if Dat.optParam.showits == 1
                fprintf(Dat.VAL.PF, Dat.VAL.itctr, Dat.VAL.Fmin, toc, Dat.VAL.nLocSearch, Dat.VAL.feval + Dat.VAL.fevalLocal);
            end
            Item_History(end + 1, :) = [Item_Iter, Dat.VAL.feval + Dat.VAL.fevalLocal, Dat.VAL.Fmin, toc];
        end

        function Stop(Dat)
            if Dat.VAL.Fmin < Dat.optParam.goal
                fprintf('The algorithm found a better solution than the target value given: %10.05f', Dat.VAL.Fmin);
                disp(' ')
            end
            if Dat.VAL.itctr >= Dat.optParam.maxits
                fprintf('A maximum number of iterations has been exceeded by the algorithm: %4i', Dat.optParam.maxits);
                disp(' ')
            end
            if Dat.VAL.feval + Dat.VAL.fevalLocal > Dat.optParam.maxevals 
                fprintf('A maximum number of function evaluations has been exceeded by the algorithm: %8i\n', Dat.optParam.maxevals);
                disp(' ')
            end
        end

        function Dat = POH_Selection( Dat )
            % Sets of Global/Local POHs
            [ Dat.POH, Dat.POHLocal ] = deal( [] );

            % Find Global set of POHs
            [ set, Dat ] = Identify_set( Dat );
            [ Dat.POH, Dat.POHLocal ] = deal( set(1, :) );
        end

    end
    
end