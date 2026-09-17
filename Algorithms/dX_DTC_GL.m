function [Fmin, Xmin, History] = dX_DTC_GL(Problem, opts, Bounds)
%--------------------------------------------------------------------------
% This function implements a deterministic global optimization algorithm 
% of the DIRECT-type. An algorithm designed to solve box-constrained continuous
% optimization problems without utilizing derivatives.
%
% Objective:
%   minimize    fun(x)
%       s.t.    lb <= x <= ub,
%
% INPUTS:
%   Problem.f - Function handle containing the objective function to minimize
%   opts      - Structure with optional user-defined settings
%               + opts.maxevals: Maximum number of function evaluations;
%                   default: 1e4
%               + opts.maxits: Maximum number of algorithm iterations; 
%                   default: 1e4
%               + opts.tolx: Termination tolerance on subdivision of hyper-rectangles 
%                   depth, a positive scalar; default: 1e-16
%               + opts.showits: Display iteration status; default: 1
%               + opts.dimension: Problem dimensionality; default: 1
%               + opts.globalmin: Known global minimum; default: -inf
%               + opts.ftarget: Target function value; default: -inf              
%               + opts.testflag: Flag for testing, 0 off, 1 looks for opts.tol 
%                   error for the given opts.globalmin, 2 aims to find better 
%                   value tha opts.ftarget; default: 0
%               + opts.tol: tolerance for termination if opts.testflag = 1;
%                   default: 1e-2
%               + opts.model: linear of quadratic approximation used to investigate 
%                   the potentiality within the hyper-rectangles; default: quadratic
%   Bounds    - Matrix specifying the bounds for the problem variables [lb, ub], 
%               where lb and ub are column vectors of the same length as 'x', 
%               representing lower and upper bounds in the constraint lb <= x <= ub.
%
% OUTPUTS:
%   Fmin    - Minimum function value found
%   Xmin    - The input value that gives Fmin
%   History - Optional output containing history of iterations
%
% REFERENES:
%   L. Stripinis, R. Paulavičius. (2026). A practical DIRECT-type algorithm for 
%   medium-scale black box global optimization.Applied Soft Computing. ISSN: 1568-4946,
%   Online first, article 116392, 46 pages. DOI: 10.1016/j.asoc.2026.116392
%--------------------------------------------------------------------------

% Handle missing input arguments
if nargin == 2, Bounds = []; end             % If only Problem and opts are given
if nargin == 1, Bounds = []; opts = []; end  % If only Problem is given

% Get options and initialize values
VAL = Options(opts, nargout, Problem, Bounds);

% Allocate sets and initialize variables
[MSS, VAL] = Initialization(VAL, Problem);

%Main optimization loop
while VAL.exitflag > 0
    % Selection of potential optimal hyper-rectangles
    [Main, POH] = Selection(VAL);

    % Subdivide potential optimal hyper-rectangles
    [VAL, MSS] = Subdivision(VAL, MSS, POH, Main);

    % Update minima and check stopping conditions
    VAL = Arewedone(VAL);
end

% Return results
Fmin = VAL.Fmin;  % Return the minimum value found
Xmin = (VAL.xU - VAL.xL) .* VAL.Xmin + VAL.xL;  % Scale the solution back to the original bounds

if VAL.options.G_nargout == 3
    History = VAL.history;  % Return history of function evaluations
end

% Print the final results
fprintf('========================= Optimization Results =========================\n');
fprintf('Optimal solution found: \n');
fprintf('x = [ \n');
fprintf('%.10f \n', Xmin);
fprintf(']\n');
fprintf('Objective function value at optimal solution: f(x) = %.10f\n', Fmin);
fprintf('Number of iterations: %d\n', History(end, 1));
fprintf('Number of function evaluations: %d\n', History(end, 2));
fprintf('Computation time: %.4f seconds\n', History(end, 4));

switch VAL.exitflag
    case -1 
        fprintf('Exit condition: Optimization terminated successfully, minima was found with Tolerance: %4i\n', VAL.options.TOL);
    case -2
        fprintf('Exit condition: Optimization terminated successfully, better minima was found than the give function target value: %15.10f \n', VAL.Fmin);
    case -3
        fprintf('Exit condition: Exceeded max iterations. Increase maxits.\n');
    case -4
        fprintf('Exit condition: Exceeded max function evaluations. Increase maxevals.\n');
end
%--------------------------------------------------------------------------
end

function VAL = Options(opts, narg, Problem, bounds)
%--------------------------------------------------------------------------
% This function sets up and initializes options for an optimization routine.
% It retrieves values for options such as maximum iterations, tolerances, 
% function evaluation limits, and other parameters based on inputs or defaults.
%
% INPUTS:
%   opts    - Structure with optional user-specified options
%   narg    - Number of output arguments expected (G_nargout)
%   Problem - Function handle for the problem containing bounds information
%   bounds  - User-provided bounds for the optimization variables
%
% OUTPUT:
%   VAL     - Structure containing the options and bounds
%--------------------------------------------------------------------------

    % If opts is not provided or empty, initialize as an empty struct
    if nargin < 3 || isempty(opts)
        opts = [];
    end

    % Assign default values for options if not provided in opts
    getOpts(opts, ...
        'maxits',       1e4, ...           % Maximum iterations
        'maxevals',     1e4, ...           % Maximum function evaluations
        'testflag',     0, ...             % Flag for test
        'tol',          1e-2, ...          % Tolerance for minimum error
        'showits',      1, ...             % Show iteration stats
        'dimension',    1, ...             % Problem dimensionality
        'globalmin',    -inf, ...          % Known global minimum (if any)
        'tolx',         1e-16, ...         % Tolerance for variable convergence
        'ftarget',      -inf, ...          % Target function value to reach
        'model',        'quadratic'...     % Approximation in hyper-rectangles
        ); 

    % Handle the case where bounds are not provided
    if isempty(bounds)
        % Get problem information from the provided problem function
        getInfo = feval(Problem.f);

        % Determine dimensionality of the problem from getInfo
        if getInfo.nx == 0
            VAL.n = dimension;  % Use the provided dimension if nx is zero
        else
            VAL.n = getInfo.nx; % Use the dimensionality from the problem definition
        end

        % Set variable bounds from the problem information
        VAL.xL = getInfo.xl(VAL.n);  % Lower bounds
        VAL.xU = getInfo.xu(VAL.n);  % Upper bounds

        % If testflag is set, retrieve the known global minimum
        if testflag == 1
            globalmin = getInfo.fmin(VAL.n);
        end
    else
        % Use user-provided bounds
        VAL.xL = bounds(:, 1);  % Left bound
        VAL.xU = bounds(:, 2);  % Right bound
        VAL.n  = size(bounds, 1);  % Dimension of the problem
    end

    % Store options in the VAL structure
    VAL.options.G_nargout   = narg;
    VAL.options.MAXits      = maxits;
    VAL.options.MAXevals    = maxevals;
    VAL.options.TESTflag    = testflag;
    VAL.options.showITS     = showits;
    VAL.options.TOL         = tol;
    VAL.options.tolx        = tolx;
    VAL.options.targetValue = ftarget;
    VAL.options.globalMIN   = globalmin;

    if strcmp('quadratic', model)
        VAL.Model = @DynamicSubQuad;
    else
        VAL.Model = @DynamicSubLin;
    end
end

function [MSS, VAL] = Initialization(VAL, Problem)
%--------------------------------------------------------------------------
% Function to allocate memory for main sets, create initial values, and
% compute useful variables for optimization routines.
%
% INPUTS:
% VAL - Structure containing parameters and configuration options
% Problem - Structure containing the function
%
% OUTPUTS:
% MSS   - Struct containing preallocated memory for optimization sets.
% VAL   - Updated VAL structure with initialized fields.
%--------------------------------------------------------------------------

    % Start timer
    VAL.time = tic;   
    
    % Compute the number of levels for exploitation based on tolerance
    VAL.level = VAL.n * ceil(log(VAL.options.tolx / max(VAL.xU - VAL.xL)) / log(1/3));
    VAL.limit = VAL.n * ceil(log(1e-10 / max(VAL.xU - VAL.xL)) / log(1/3));
    VAL.CE = zeros(1, VAL.level + 1);

    % Initialize iteration and evaluation counters
    [VAL.evals, VAL.CE(1), VAL.exitflag, VAL.POH] = deal(1); 
    [VAL.numLsEvls, VAL.numLs, VAL.iter, VAL.succesLs] = deal(0);
    VAL.Splits = zeros(VAL.n, 1);
    [VAL.MV, VAL.MD] = deal(ones(4, 1)); 

    % Fields for the MSS structure and their corresponding sizes
    VAL.fields = {'L', 'C', 'X', 'F', 'E'};
    VAL.sizes = [VAL.n, VAL.n, VAL.n, 1, 1];
    
    % Determine functions
    VAL.fhd = @(x) feval(Problem.f, (VAL.xU - VAL.xL).*x + VAL.xL);
    VAL.Problem = @(x) feval(Problem.f, x);
    VAL.print = 'Iteration: %4i f(x) evaluations: %6i Fmin: %15.10f time(s): %10.05f \n';
    VAL.LO = optimoptions('fmincon', ...
        'MaxFunEvals',1e3*VAL.n, ...
        'MaxIter',1e3*VAL.n, ...
        'Display','off', ...
        'Algorithm','sqp', ...
        'StepTolerance',1e-10,...
        'OptimalityTolerance',1e-6,...
        'FiniteDifferenceType', 'central' ...
        );
    
    % Allocate MSS structure with preallocated fields
    MSS = struct('F', VAL.fhd(ones(VAL.n, 1)/2), ... % Preallocate F (scalar)
                 'E', 1, ...                         % Preallocate E (scalar)
                 'C', ones(VAL.n, 1)/2, ...          % Preallocate C (n x 1 vector)
                 'L', ones(VAL.n, 1)/2);             % Preallocate L (n x 1 vector)

    [VAL.Fmin, VAL.FminDirect, VAL.Median] = deal(MSS(1).F(1));
    VAL.Xmin = MSS(1).C(:, 1);
    
    % aggressive start
    LO = optimoptions('fmincon', ...
        'MaxFunEvals',1e2*VAL.n, ...
        'MaxIter',1e2*VAL.n, ...
        'Display','off', ...
        'Algorithm','sqp', ...
        'StepTolerance',1e-10,...
        'OptimalityTolerance',1e-6,...
        'FiniteDifferenceType', 'forward' ...
        );

    VAL = RunLocal( abs(VAL.xU - VAL.xL).*MSS.C(:, 1) + VAL.xL, VAL, LO );

    % Allocate history if necessary
    if VAL.options.G_nargout == 3
        VAL.history = [VAL.iter, VAL.evals, VAL.Fmin, toc(VAL.time)];
    end
    if VAL.options.showITS == 1
        fprintf(VAL.print, VAL.history);
    end
end

function [Main, POH] = Selection(VAL)
%--------------------------------------------------------------------------
% This function identifies the set of potential optimal hyper-rectangles
% (POH) based on the minimum function values (F) and distances (D).
%
% INPUTS:
%   VAL.MV - Matrix containing function values and related information for set A
%   VAL.MD - Matrix containing distances and related information for set B
%
% OUTPUTS:
%   Main - A matrix containing the indices and values of potential optimal 
%          hyper-rectangles.
%   PH   - A cell array containing the union of POH sets from F and D.
%--------------------------------------------------------------------------

    % Initialize counters and placeholders for POH sets and indices
    level = min(VAL.level, size(VAL.MV, 2));
    fmin = VAL.MV(1, 1:level);                  
    dmin = VAL.MD(1, 1:level);               
    POH = cell(1, size(fmin, 2)); 
    [index_a, index_b] = deal(size(fmin, 2)); 
    
    % Find index set of potential optimal hyper-rectangles based on F
    [m_m, index_a] = min(fmin(1:index_a));
    ss_a(1) = index_a;      
    POH{index_a} = VAL.MV(2, index_a);
    index_a = index_a - 1;
    while index_a ~= 0 && ~isnan(m_m)
        [m_m, index_a] = min(fmin(1:index_a));
        if ~isnan(m_m)                      
            ss_a(end + 1) = index_a;            %#ok<*AGROW> 
            POH{index_a} = VAL.MV(2, index_a);
        end
        index_a = index_a - 1;              
    end

    % Find index set of potential optimal hyper-rectangles based on D
    [m_m, index_b] = min(dmin(1:index_b));
    ss_b(1) = index_b;      
    POH{index_b} = union(POH{index_b}, VAL.MD(2, index_b));
    index_b = index_b - 1;
    while index_b ~= 0 && ~isnan(m_m)
        [m_m, index_b] = min(dmin(1:index_b));
        if ~isnan(m_m)                     
            ss_b(end + 1) = index_b;            
            POH{index_b} = union(POH{index_b}, VAL.MD(2, index_b)); 
        end
        index_b = index_b - 1;            
    end

    % Remove indices from ss_b that are also in the intersection of MM(3, ss_a) and DD(3, ss_b)
    common_indices = intersect(VAL.MD(3, ss_b), VAL.MV(3, ss_a));
    ss_b(ismember(VAL.MD(3, ss_b), common_indices)) = [];

    % Create the Main matrix
    Main = [ss_a, ss_b; VAL.MV(2, ss_a), VAL.MD(2, ss_b); VAL.MV(3, ss_a), VAL.MD(3, ss_b); VAL.MV(4, ss_a), VAL.MD(4, ss_b); zeros(1, length(ss_a)), ones(1, length(ss_b))];
    [~, sortOrder] = sort(Main(3, :), 'descend');
    Main = Main(:, sortOrder);
    [~, sortOrder] = sort(Main(1, :), 'descend');
    Main = Main(:, sortOrder);
end

function [VAL, MSS] = Subdivision(VAL, MSS, POH, Main)
%--------------------------------------------------------------------------
% This function manages the storage of MSS data, processes the POH (potential 
% optimal hyper-rectangles), and cleans up data by removing entries from MSS 
% according to the POH. After processing, it updates the minimum values in VAL.
%
% INPUTS:
%   VAL  - Structure containing evaluation counters and minimum values
%   MSS  - Structure containing fields 'F', 'E', 'L', 'C' (main storage sets)
%   POH  - Cell array of potential optimal hyper-rectangles for each MSS index
%   Main - Data matrix for partitioning
%
% OUTPUTS:
%   VAL  - Updated structure with new evaluation data and minimum values
%   MSS  - Updated storage structure with processed entries
%--------------------------------------------------------------------------

    % Update MSS with the new data from Main using the Storage function
    [MSS, VAL, POH] = Storage(MSS, VAL, Main, POH);

    % Loop backwards through POH 
    for i = size(POH, 2):-1:1
        if ~isempty(POH{i}) 
            
            % Check if the size of POH{i} matches the number of evaluations in VAL.CE(i)
            if (VAL.CE(i) - size(POH{i}, 2)) == 0
                if find(VAL.CE ~= 0, 1, 'first') == i
                    [MSS(i).E, MSS(i).L, MSS(i).F, MSS(i).C, MSS(i).X] = deal([]);
                end
            else
                C = 1:VAL.CE(i);  % Indices of current evaluations
                C(POH{i}) = [];   % Exclude the indices in POH{i}
                
                % Prepare the new positions for the remaining elements
                pp = min(POH{i}):length(C);
                
                % Reassign the non-POH elements to their new positions in MSS(i)
                MSS(i).E(pp) = MSS(i).E(C(pp));
                MSS(i).L(:, pp) = MSS(i).L(:, C(pp));
                MSS(i).F(pp) = MSS(i).F(C(pp));
                MSS(i).C(:, pp) = MSS(i).C(:, C(pp));
            end
            
            % Update the count of evaluations after removal
            VAL.CE(i) = VAL.CE(i) - size(POH{i}, 2);
        end
    end

    % After processing POH, update the minimum values in VAL
    VAL = Find_min(MSS, VAL);
end

function [MSS, VAL, POH] = Storage(MSS, VAL, Main, POH)
%--------------------------------------------------------------------------
% This function updates the MSS structure with new partition data from Main.
% It checks if new storage allocation is required for each field in MSS 
% and updates fields based on the calculated partition (DD) from the Partitioning function.
%
% INPUTS:
%   MSS   - Structure containing fields 'F', 'E', 'L', 'C' (main storage sets)
%   VAL   - Structure containing evaluation counters, sizes, fields, and CE (storage tracking)
%   Main  - Matrix that provides indices for partitioning and updating MSS
%
% OUTPUTS:
%   MSS   - Updated MSS structure with the newly added partition data
%   VAL   - Updated VAL structure with updated counters and storage information
%--------------------------------------------------------------------------

    % Loop through each column of Main
    % fminOld = VAL.Xmin;
    for i = 1:size(Main, 2)
        Rule = true;
        VAL.frox = inf(VAL.n, 1);
        while Rule
            % Partition the data using the Partitioning function
            [DD, mdx, VAL, ls] = Partitioning(i, MSS, VAL, Main);
    
            % Check if new memory allocation is needed for MSS(mdx)
            if VAL.CE(mdx) == 0
                % Initialize storage for each field
                for j = 1:length(VAL.fields)
                    MSS(mdx).(VAL.fields{j}) = zeros(VAL.sizes(j), VAL.n * 10);
                end
            elseif VAL.CE(mdx) > size(MSS(mdx).F, 2) + 3
                % Expand storage if current storage is insufficient
                for j = 1:length(VAL.fields)
                    MSS(mdx).(VAL.fields{j}) = [MSS(mdx).(VAL.fields{j}), zeros(VAL.sizes(j), VAL.n * 10)];
                end
            end
    
            % Assign indices for new entries (3 new entries)
            II = [VAL.CE(mdx) + 1, VAL.CE(mdx) + 2, VAL.CE(mdx) + 3];
            VAL.CE(mdx) = II(3);  % Update the storage counter
    
            % Update MSS fields with the new partition data (DD)
            MSS(mdx).F(II) = DD.F;
            MSS(mdx).E(II) = DD.E;
            MSS(mdx).L(:, II) = DD.L .* ones(1, 3);
            MSS(mdx).C(:, II) = DD.C;
            
            % Update Fmin
            VAL.Fmin = min([DD.F, VAL.Fmin]);
            
            [id, VAL, ImRul] = VAL.Model(DD.F, DD.C, DD.E, VAL, ls);
            VAL.FminDirect = min([DD.F, VAL.FminDirect]);

            Rule = ImRul && VAL.limit > mdx;
            % Rule = false;
            if Rule
                if mdx > numel(POH)
                    POH{mdx} = II(id);
                else
                    POH{mdx}(end + 1) = II(id); 
                end
            end
            Main(:, i) = [mdx; II(id); DD.E(id); DD.F(id); Main(5, i)];
        end
    end
end

function [A, mdx, VAL, ls] = Partitioning(i, MSS, VAL, Main)
%--------------------------------------------------------------------------
% This function calculates and updates certain properties (A.C, A.F, A.L, A.E)
% based on the input values from MSS, VAL, and Main. It updates the working 
% set and increments evaluations accordingly.
%
% INPUTS:
%   i     - The index used to select the main sets from MSS and Main
%   MSS   - A structure containing sets 'C', 'F', 'L', 'E'
%   VAL   - A structure with necessary parameters 
%   Main  - A matrix for indexing into MSS (contains sets and index values)
%
% OUTPUTS:
%   A     - A structure containing updated values for 'C', 'F', 'L', 'E'
%   mdx   - Updated index after processing
%   VAL   - Updated VAL structure
%--------------------------------------------------------------------------

    % Retrieve the main set index for MSS based on Main
    mdx = Main(1, i);
    
    % Extract and replicate values from MSS into A structure
    A.C = MSS(mdx).C(:, Main(2, i)) .* ones(1, 3);  % Replicate column 3 times
    A.F = MSS(mdx).F(Main(2, i)) .* ones(1, 3);     % Replicate scalar 3 times
    A.L = MSS(mdx).L(:, Main(2, i));                % Take as-is (n x 1 vector)
    A.E = MSS(mdx).E(Main(2, i)) .* ones(1, 3);     % Replicate scalar 3 times

    % Find the index where A.L is maximized
    ls = find(A.L == max(A.L));
    delta = A.L(ls(1)) * (2/3);
    ls = ls(find(VAL.Splits(ls) == min(VAL.Splits(ls)), 1, "first"));
    VAL.Splits(ls) = VAL.Splits(ls) + 1;

    % Update A.C and A.F for the new DELTA-adjusted values
    A.C(ls, 2) = A.C(ls, 2) - delta;  % Adjust C for the first value
    A.F(2) = VAL.fhd(A.C(:, 2));

    A.C(ls, 3) = A.C(ls, 3) + delta;  % Adjust C for the second value
    A.F(3) = VAL.fhd(A.C(:, 3));

    % Update evaluation counters and A.E
    A.E(2) = VAL.evals + 1; 
    A.E(3) = VAL.evals + 2; 
    VAL.evals = A.E(3);  % Increment VAL.evals by 2

    % Update A.L at the selected index
    A.L(ls) = delta / 2;

    % Increment mdx to reflect the updated set index
    mdx = mdx + 1;
    %--------------------------------------------------------------------------
end

function VAL = Arewedone(VAL)
%--------------------------------------------------------------------------
% This function checks the termination conditions of an optimization routine.
% It checks if the minimum has been found within a specified tolerance or 
% if any iteration or evaluation limits have been exceeded. It also stores
% iteration history and optionally displays iteration statistics.
%
% INPUT:
%   VAL - Structure containing the current state of the optimization
%
% OUTPUT:
%   VAL - Updated structure with potential termination flag or updated history
%--------------------------------------------------------------------------
    % Increment the iteration number
    VAL.iter = VAL.iter + 1;

    % Store iteration history if required (G_nargout == 3 indicates history storage)
    if VAL.options.G_nargout == 3
        % Append current iteration stats to history (iter, evals, Fmin, elapsed time)
        VAL.history(end + 1, :) = [VAL.iter, VAL.evals + VAL.numLsEvls, VAL.Fmin, toc(VAL.time)];
    end

    % Optionally display iteration statistics
    if VAL.options.showITS == 1
        % Print the last entry in history (iteration stats)
        fprintf(VAL.print, VAL.history(end, :));
    end

    % Check if a known global minimum is provided and calculate error
    if VAL.options.TESTflag == 1
        % Calculate the percentage error if a global minimum is known
        if VAL.options.globalMIN ~= 0
            VAL.exitflag = 100 * (VAL.Fmin - VAL.options.globalMIN) / abs(VAL.options.globalMIN);
        else
            VAL.exitflag = 100 * VAL.Fmin;
        end

        % Check if the error is within the specified tolerance
        if VAL.exitflag < VAL.options.TOL
            VAL.exitflag = -1;
        end
    elseif VAL.options.TESTflag == 2
        if VAL.options.targetValue >= VAL.Fmin
            VAL.exitflag = -2;
        end
    end

    % Check if the maximum number of iterations has been reached
    if VAL.iter >= VAL.options.MAXits
        VAL.exitflag = -3;
    end

    % Check if the maximum number of function evaluations has been reached
    if VAL.evals + VAL.numLsEvls > VAL.options.MAXevals
        VAL.exitflag = -4;
    end
end

function VAL = Find_min(MSS, VAL)
%--------------------------------------------------------------------------
% This function finds the minimum function value (Fmin) and its corresponding
% location (Xmin) from the MSS structure. Additionally, it computes the
% Euclidean distance between Xmin and other points in MSS.
%
% INPUTS:
%   MSS - Structure containing fields F (function values), E (evaluations), and C (coordinates)
%   VAL - Structure containing current evaluation counters and will store min values and distances
%
% OUTPUT:
%   VAL - Updated structure containing:
%         MV  - Minimum values (function values, indices, evaluations)
%         MD  - Distances (Euclidean distance, indices, evaluations)
%         Fmin, Xmin - The minimum function value and corresponding coordinates
%--------------------------------------------------------------------------
    
    % Initialize VAL.MV and VAL.MD with NaN values for storing results
    [VAL.MV, VAL.MD] = deal(nan(4, size(MSS, 2)));

    % Loop over each column of MSS to find minimum function values
    for i = 1:size(MSS, 2)
        if VAL.CE(i) ~= 0
            % Find the index of the maximum evaluation counter corresponding to the minimum function value (F)
            minF = min(MSS(i).F(1:VAL.CE(i)));  % Find minimum function value in current column
            maxEvalIdx = find(MSS(i).E(1:VAL.CE(i)) == max(MSS(i).E(MSS(i).F(1:VAL.CE(i)) == minF)));

            % Store the minimum function value, its index, and the corresponding evaluation
            VAL.MV(1, i) = MSS(i).F(maxEvalIdx);  % Function value
            VAL.MV(2, i) = maxEvalIdx;            % Index
            VAL.MV(3, i) = MSS(i).E(maxEvalIdx);  % Evaluation
            VAL.MV(4, i) = MSS(i).F(maxEvalIdx);  % Function value
        end
    end
    
    % Find the overall minimum function value (Fmin) and its index
    [VAL.FminDirect, fminindex] = min(VAL.MV(1, :));
    
    if VAL.FminDirect <= VAL.Fmin
        VAL.Fmin = VAL.FminDirect;
        VAL.Xmin = MSS(fminindex).C(:, VAL.MV(2, fminindex));
    end

    % Calculate the Euclidean distance between Xmin and all points in MSS
    for i = 1:size(MSS, 2)
        if VAL.CE(i) ~= 0
            % Calculate Euclidean distances between Xmin and all points in column
            D = sqrt(sum((VAL.Xmin - MSS(i).C(:, 1:VAL.CE(i))).^2, 1));

            % Find the index of the maximum evaluation counter corresponding to the minimum distance
            minDistIdx = find(MSS(i).E(1:VAL.CE(i)) == max(MSS(i).E(D == min(D))));

            % Store the minimum distance, its index, and the corresponding evaluation
            VAL.MD(1, i) = D(minDistIdx);         % Distance
            VAL.MD(2, i) = minDistIdx;            % Index
            VAL.MD(3, i) = MSS(i).E(minDistIdx);  % Evaluation
            VAL.MD(4, i) = MSS(i).F(minDistIdx);  % Evaluation
        end
    end
end

function [POH, VAL, Rule] = DynamicSubQuad(F, C, E, VAL, ls)
%--------------------------------------------------------------------------
% Function: DynamicSub
% Purpose: Determine if the hyper-rectangle will be subdivided further
% Inputs:
%   datStr: current algorithm samples
%   index: indices of points
%   opts: structure with options and status of optimization processes
% Outputs:
%   POH: target hyper-rectangle
%   opts: updated structure
%   Rule: decision on hyper-rectangle
%--------------------------------------------------------------------------

    % Calculate quadratic coefficients
    G = 3;
    X = [0; -2; 2]; 
    X_vals = [ones( 3, 1 ), X, X.^2];
    b_vals = X_vals \ F';

    % Determine optimum
    if b_vals(3) ~= 0
        xopt = max(-G, min(G, -b_vals(2)/(2*b_vals(3))));
    else
        xopt = 0;
    end
    points = [-G; xopt; G]; 
    y0 = [ones(3,1), points, points.^2]*b_vals; 
    [f_min, f_id] = min(y0);

    % Better than Fmin?
    [~, POH] = min( abs( points(f_id) - X ) );
    Rule = f_min - VAL.Fmin + 1e-14 < 0;

    % Bound check
    delta = C(ls, 3) - C(ls, 1);
    if round(C(ls, 3) + delta, 16) == 1 && round(C(ls, 2) - delta, 16) == 0
        Bounds = points(f_id) <= G && points(f_id) >= -G;
    elseif round(C(ls, 2) - delta, 16) == 0
        Bounds = points(f_id) < G && points(f_id) >= -G;
    elseif round(C(ls, 3) + delta, 16) == 1
        Bounds = points(f_id) <= G && points(f_id) > -G;
    else
        Bounds = points(f_id) < G && points(f_id) > -G;
    end

    if ~ismember(E(POH), VAL.POH) && ((Rule && Bounds) || VAL.FminDirect - 1e-8 > min(F))
        % Starting point for local solver
        VAL.POH(end+1) = E(POH);

        % Run local solver
        VAL = RunLocal(abs(VAL.xU - VAL.xL).*C(:, POH) + VAL.xL, VAL, VAL.LO);
    end
end

function [POH, VAL, Rule] = DynamicSubLin(F, C, E, VAL, ls)
%--------------------------------------------------------------------------
% Function: DynamicSub
% Purpose: Determine if the hyper-rectangle will be subdivided further
% Inputs:
%   datStr: current algorithm samples
%   index: indices of points
%   opts: structure with options and status of optimization processes
%   ls: indice of subdivided coordinate 
%   delta: distance between points
% Outputs:
%   POH: target hyper-rectangle
%   opts: updated structure
%   Rule: decision on hyper-rectangle
%--------------------------------------------------------------------------

    % Design left linear model
    X1 = [ones(2, 1), [0; -2]]; 
    b1 = ((X1'*X1)\X1')*F([1,2])'; 
    MetaF1 = @(x) ([1, x])*b1;
    
    % Design right linear model
    X2 = [ones(2, 1), [0; 2]];  
    b2 = ((X2'*X2)\X2')*F([1,3])'; 
    MetaF2 = @(x) ([1, x]) * b2;
    
    fm = [MetaF2(-1), MetaF1(1), MetaF1(-3), MetaF2(3)]; 
    points = [0, -3, 3];

    [f_min, POH] = min([min(fm(1:2)), fm(3), fm(4)]); 
    Rule = f_min - VAL.Fmin + 1e-16 < 0;

    % Bound check
    delta = C(ls, 3) - C(ls, 1);
    if round(C(ls, 3) + delta, 16) == 1 && round(C(ls, 2) - delta, 16) == 0
        Bounds = points(POH) <= 3 && points(POH) >= -3;
    elseif round(C(ls, 2) - delta, 16) == 0
        Bounds = points(POH) < 3 && points(POH) >= -3;
    elseif round(C(ls, 3) + delta, 16) == 1
        Bounds = points(POH) <= 3 && points(POH) > -3;
    else
        Bounds = points(POH) < 3 && points(POH) > -3;
    end

    if ~ismember(E(POH), VAL.POH) && ((Rule && Bounds) || VAL.FminDirect - 1e-8 > min(F))
        VAL.POH(end+1) = E(POH);

        % Run local solver
        VAL = RunLocal( abs(VAL.xU - VAL.xL).*C(:, POH) + VAL.xL, VAL, VAL.LO );
    end
end

function VAL = RunLocal(x0, VAL, opts)
%--------------------------------------------------------------------------
% Function: RunLocal
% Purpose: Runs local solver from any given point
% Inputs:
%   opts: structure with options and status of optimization processes
%   x0: starting point for local solver
% Outputs:
%   opts: updated structure
%--------------------------------------------------------------------------
% Run Local Solver
    [x, fval, ~, output] = fmincon(VAL.Problem, x0, [], [], [], [], VAL.xL,...
        VAL.xU, [], opts);

    % Update optimization stats
    VAL.numLs = VAL.numLs + 1;
    VAL.numLsEvls = VAL.numLsEvls + output.funcCount;
    
    % Update global minimum if improved
    if fval < VAL.Fmin
        VAL.succesLs = VAL.succesLs + 1;
        VAL.Fmin = fval;
        VAL.Xmin = (x - VAL.xL)./(VAL.xU - VAL.xL);
    end
end

function varargout = getOpts(options, varargin)
%--------------------------------------------------------------------------
% Function: getOpts
% Purpose: Returns options values in an options structure
% Inputs:
%   options: structure with options for optimization
%   varargin: default values for optimization
% Outputs:
%   varargout: parameters for optimization
%--------------------------------------------------------------------------
    K = fix(nargin/2);
    varargout = cell(K, 1);
    for i = 1:K
        if isa(options, 'struct') && isfield(options, varargin{i*2 - 1})
            assignin('caller', varargin{i*2 - 1}, options.(varargin{i*2 - 1}));
        else
            assignin('caller', varargin{i*2 - 1}, varargin{i * 2});
        end
    end
end
%--------------------------------------------------------------------------
% END of BLOCK
%--------------------------------------------------------------------------