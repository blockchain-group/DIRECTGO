function [minima, xatmin, history] = dDirect_NAS(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : dDirect_NAS
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dDirect_NAS(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       NAS       - Neighbourhood Assignment Strategy for
%                   problems with constraints
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f           = Objective function handle
%
%       bounds  - (n x 2) matrix of bound constraints LB <= x <= UB
%                 The first column is the LB bounds and the second
%                 column contains the UB bounds
%
%       opts    - MATLAB structure which contains options.
%                 opts.maxevals  = max. number of function evals
%                 opts.maxits    = max. number of iterations
%                 opts.maxdeep   = max. number of rect. divisions
%                 opts.testflag  = 1 if globalmin known, 0 otherwise
%                 opts.globalmin = globalmin (if known)
%                 opts.ept       = tollerance for constrains
%                 opts.ep        = global/local weight parameter
%                 opts.showits   = 1 print iteration status
%                 opts.pert      = pertubation for implicit constraint
%                                  handling
%                 opts.tol       = tolerance for termination if
%                                  testflag = 1
%
% Output parameters:
%       minima  -  best minimum value which was founded
%
%       xatmin  - coordinate of minimal value
%
%       history - (iterations x 4) matrix of iteration history
%                 First column coresponds iterations
%                 Second column coresponds number of objective function
%                 evaluations
%                 Third column coresponds minima value of objecgtive
%                 function which was founded at iterations
%                 Third column coresponds time cost of the algorithm
%
% Original DIRECT implementation taken from:
%--------------------------------------------------------------------------
% D.R. Jones, C.D. Perttunen, and B.E. Stuckman. "Lipschitzian
% Optimization Without the Lipschitz Constant". Journal of Optimization
% Theory and Application, 79(1):157-181, (1993). DOI 10.1007/BF00941892
%
% Constraint handling strategie taken from:
%--------------------------------------------------------------------------
% J. M. Gablonsky. "Modi?cations of the DIRECT algorithm". Ph.D. thesis,
% North Carolina State University (2001).
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[SS, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[MSS, CE, third, VAL] = Alocate(SS, VAL);

% Initialization step
[MV, MSS, CE, VAL, minval, xatmin] = Initialization(Problem,...
    MSS, CE, VAL, SS);

while VAL.perror > SS.TOL                 % Main loop
%--------------------------------------------------------------------------
    % Selection of potential optimal hyper-rectangles step
    POH = Find_po(MSS, CE, SS, minval, MV(1, :));
    
    % Subdivide potential optimalhyper-rectangles
    [VAL, MSS, CE, MV, minval, xatmin] = Calulcs(VAL, Problem,...
        MSS, CE, third, POH, SS);
    
    % Update minima and check stopping conditions
    [VAL, SS] = Arewedone(minval, VAL, SS, MSS);
    
    % Replace infeasible hyper-rectangles with feasible values from the
    % neighbourhood
    MSS = replaceinf(MSS, CE, VAL, SS.pert);
    
%--------------------------------------------------------------------------
end                                       % End of while

% Return value
minima      = minval;
if SS.G_nargout == 2
    xatmin    = abs(VAL.b - VAL.a).*xatmin(:, 1) + VAL.a;
elseif SS.G_nargout == 3
    xatmin    = abs(VAL.b - VAL.a).*xatmin(:, 1) + VAL.a;
    history   = VAL.history(1:(VAL.itctr), 1:4);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% AUXILIARY FUNCTION BLOCK
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Function  : SS
% Purpose   : Get options from inputs
%--------------------------------------------------------------------------
function [OPTI, VAL, Problem] = Options(opts, narg, Problem, bounds)
%--------------------------------------------------------------------------
% Determine option values
if nargin < 3 && isempty(opts)
    opts = [];
end
getopts(opts, 'maxits', 1000,'maxevals', 100000, 'maxdeep', 1000,...
    'testflag', 0, 'tol', 0.01, 'showits', 1, 'dimension', 1, 'ept',...
    1e-8, 'globalmin', 0, 'globalxmin', 0, 'ep', 1e-4, 'pert', 1e-6);

if isempty(bounds)

% Return the problem information.
getInfo = feval(Problem.f);

if isfield(getInfo, 'confun')
    Problem.constraint = getInfo.confun;
end
    
% dimension
    if getInfo.nx == 0
        VAL.n = dimension;
    else
        VAL.n = getInfo.nx;
    end
    
    VAL.a = arrayfun(@(i) getInfo.xl(i), 1:VAL.n)'; % left bound
    VAL.b = arrayfun(@(i) getInfo.xu(i), 1:VAL.n)'; % right bound
    
    if testflag == 1
% minimum value of function
        OPTI.globalMIN  = getInfo.fmin(VAL.n);
% minimum point of function
        OPTI.globalXMIN = arrayfun(@(i) getInfo.xmin(i), 1:VAL.n)';
    end
else
    VAL.a = bounds(:, 1);               % left bound
    VAL.b = bounds(:, 2);               % right bound
    VAL.n = size(bounds, 1);            % dimension
    if testflag == 1
        OPTI.globalMIN = globalmin;     % minimum value of function
        OPTI.globalXMIN = globalxmin;   % minimum point of function
    end
end

OPTI.G_nargout = narg;     % output arguments
OPTI.MAXits    = maxits;   % Fmax of iterations
OPTI.MAXevals  = maxevals; % Fmax # of function evaluations
OPTI.MAXdeep   = maxdeep;  % Fmax number of side divisions
OPTI.TESTflag  = testflag; % terminate if global minima is known
OPTI.showITS   = showits;  % print iteration stat
OPTI.TOL       = tol;      % allowable relative error if f_reach is set
OPTI.ept       = ept;      % tollerance for constraints
OPTI.ep        = ep;       % global/local weight parameter
OPTI.pert      = pert;     % pertubation for implicit constraint handling
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, third, VAL] = Alocate(SS, VAL)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
VAL.time       = 0;                         % initial time
VAL.fcount     = 1;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
CE             = zeros(1, 10*VAL.n*SS.MAXdeep);% collum counter

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
    'L', zeros(VAL.n, 1));

third           = zeros(1, 10*SS.MAXdeep);  % delta values
third(1)        = 1/3;                      % first delta
for i = 2:10*VAL.n*SS.MAXdeep               % all delta
    third(i)    = (1/3)*third(i - 1);
end
if SS.G_nargout == 3
    VAL.history = zeros(SS.MAXits, 4);       % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [MV, MSS, CE, VAL, minval, xatmin] = Initialization(Problem,...
    MSS, CE, VAL, SS)
%--------------------------------------------------------------------------
% Create Initial values
MV(1:3, 1) = 1;                                      % first fake value
MSS(1).L(:, 1) = zeros(VAL.n, 1);                    % Lengths
MSS(1).C(:, 1) = ones(VAL.n, 1)/2;                   % Center point
MSS(1).E(1) = 1;                                     % Index
MSS(1).Diam  = 1/2*norm((1/3*(ones(VAL.n, 1))).^MSS(1).L(:, 1));
point = abs(VAL.b - VAL.a).*MSS(1).C(:,1) + VAL.a;  % Real point
if isfield(Problem, 'constraint')
    MSS.X(1) = CallC(Problem, point, SS);
    if MSS.X(1) == 0
        [MSS(1).F(1), minval, MV(1)] = deal(feval(Problem.f, point));
    else
        [MSS(1).F(1), minval, MV(1)] = deal(10^9);
    end
else
    ff = feval(Problem.f, point);
    if ~isnan(ff) && ~isinf(ff) 
        [MSS(1).F(1), minval, MV(1)] = deal(ff);
        MSS.X(1) = 0;
    else
        MSS.X(1) = 1;
        [MSS(1).F(1), minval, MV(1)] = deal(10^9);
    end
end
CE(1) = 1;                                          % \# in column.
xatmin = MSS(1).C(:, 1);

% Check stop condition if global minima is known
if SS.TESTflag  == 1
    if SS.globalMIN ~= 0
        VAL.perror = 100*(minval - SS.globalMIN)/abs(SS.globalMIN);
    else
        VAL.perror = 100*minval;
    end
else
    VAL.perror   = 2;
end

% Store History
if SS.G_nargout == 3
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.fcount;
    VAL.history(VAL.itctr, 3) = minval;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
%--------------------------------------------------------------------------
function [VAL, SS] = Arewedone(minval, VAL, SS, MSS)
%--------------------------------------------------------------------------
VAL.time = toc;

% Show iteration stats
if SS.showITS == 1
    fprintf(...
        'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, minval, VAL.time, VAL.fcount);
end

% Check for stop condition
if SS.TESTflag == 1
    % Calculate error if globalmin known
    if SS.globalMIN ~= 0
        VAL.perror = 100*(minval - SS.globalMIN)/abs(SS.globalMIN);
    else
        VAL.perror = 100*minval;
    end
    if VAL.perror < SS.TOL
        fprintf('Minima was found with Tolerance: %4i', SS.TOL);
        VAL.perror = -1;
    end
else
    VAL.perror = 10;
end

% Have we exceeded the maxits?
if VAL.itctr >= SS.MAXits
    disp('Exceeded max iterations. Increase maxits'); VAL.perror = -1;
end

% Have we exceeded the maxevals?
if VAL.fcount > SS.MAXevals
    disp('Exceeded max fcn evals. Increase maxevals'); VAL.perror = -1;
end

% Have we exceeded max deep?
if SS.MAXdeep <= max(MSS(end).L(:, 1)) + 1
    disp('Exceeded Max depth. Increse maxdeep'); VAL.perror = -1;
end

% Store History
if SS.G_nargout == 3
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.fcount;
    VAL.history(VAL.itctr, 3) = minval;
    VAL.history(VAL.itctr, 4) = VAL.time;
end

% Update iteration number
if VAL.perror > SS.TOL
    VAL.itctr = VAL.itctr + 1;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [VAL, MSS, CE, MV, minval, xatmin] = Calulcs(VAL,...
    Problem, MSS, CE, third, POH, OPTI)
%--------------------------------------------------------------------------
[MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, VAL, POH, OPTI);

for i = size(POH, 2):-1:1                
    if ~isempty(POH{i})         
        if (CE(i) - size(POH{i}, 2)) == 0
            if find(CE ~= 0, 1, 'first') == i
                [MSS(i).E, MSS(i).L, MSS(i).F, MSS(i).C, MSS(i).X] = deal([]);
            end
        else
            C = setdiff(1:CE(i), POH{i});
            MSS(i).E(1:length(C)) = MSS(i).E(C);
            MSS(i).X(1:length(C)) = MSS(i).X(C);
            MSS(i).L(:, 1:length(C)) = MSS(i).L(:, C);
            MSS(i).F(1:length(C)) = MSS(i).F(C);
            MSS(i).C(:, 1:length(C)) = MSS(i).C(:, C);
        end
        CE(i) = CE(i) - size(POH{i}, 2);
    end
end

% Find minima values
[MV, minval, xatmin] = Find_min(CE, MSS);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, VAL, PH, OPTI)
%--------------------------------------------------------------------------
for i = 1:size(MSS, 2)
    if ~isempty(PH{i})
        for j = 1:size(PH{i}, 2)
            [DD, TT, VAL] = Evaluations(Problem, i, MSS, third,...
                VAL, PH{i}(j), OPTI);
            for h = 1:TT
                II              = i + h;               % SET index
                if II > length(CE)
                    CE(II) = 0;
                end
                if CE(II) == 0
                    MSS(II).Diam = 1/2*norm((1/3*(ones(VAL.n,...
                                    1))).^DD.lL(:, h));
                    [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n,...
                                                CE(II - 1)*2));
                    [MSS(II).F, MSS(II).E, MSS(II).X] = deal(zeros(1,...
                                                        CE(II - 1)*2));
                    MSS(II).X = nan(1, 100);
                end
                IL              = CE(II) + 1;            % Left index
                IR              = CE(II) + 2;            % Right index
                CE(II)          = IR;                    % Colum index
                if CE(II) > size(MSS(II).F, 2)
                    MSS(II).L = [MSS(II).L, zeros(VAL.n, CE(II - 1)*2)];
                    MSS(II).F  = [MSS(II).F, nan(1, CE(II - 1)*2)];
                    MSS(II).E = [MSS(II).E, zeros(1, CE(II - 1)*2)];
                    MSS(II).C = [MSS(II).C, zeros(VAL.n, CE(II - 1)*2)];
                    MSS(II).X = [MSS(II).X, nan(1, CE(II - 1)*2)];
                end
                MSS(II).F(IL)    = DD.L(DD.lsx(h, 1));   % Left f(x)
                MSS(II).F(IR)    = DD.R(DD.lsx(h, 1));   % Right f(x)
                MSS(II).X(IL)    = DD.Lx(DD.lsx(h, 1));  % Left flag
                MSS(II).X(IR)    = DD.Rx(DD.lsx(h, 1));  % Right flag
                MSS(II).E(IL)    = DD.eL(h);             % Left fcn counter
                MSS(II).E(IR)    = DD.eR(h);             % Right fcn counter
                MSS(II).L(:, IL) = DD.lL(:, h);          % Left lenght
                MSS(II).L(:, IR) = DD.lL(:, h);          % Right lenght
                MSS(II).C(:, IL) = DD.cL(:, DD.lsx(h, 1)); % Left x
                MSS(II).C(:, IR) = DD.cR(:, DD.lsx(h, 1)); % Right x
            end
            II                    = i + TT;            % SET index
            CE(II)                = CE(II) + 1;        % Colum index
            MSS(II).F(CE(II))     = DD.O;              % Center f(x)
            MSS(II).X(CE(II))     = DD.Ox;             % Center flag
            MSS(II).E(CE(II))     = DD.eO;             % Center fcn counter
            MSS(II).C(:, CE(II))  = DD.cO;            % Center x
            MSS(II).L(:, CE(II))  = DD.LO;            % Center lenght
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function final_pos = Find_po(MSS, CE, SS, minval, MV)
% Find all rects on hub

index    = size(MSS, 2);
nonsize  = [MSS.Diam];
hulls    = cell(3, index);

for i = 1:index
    if CE(i) ~= 0
        ties     = find(abs(MSS(i).F(1:CE(i)) - MV(i)) <= 1e-13);
        if ~isempty(ties)
            hulls{1, i} = ties;
            hulls{2, i} = nonsize(i)*ones(1, size(ties, 2));
            hulls{3, i} = MSS(i).F(ties);
        end
    end
end

fc   = cell2mat(hulls(3, :));
szes = cell2mat(hulls(2, :));

% Compute lb and ub for rects on hub
ubound = calc_ubound(fc, hulls, szes, index);
lbound = calc_lbound(fc, hulls, szes, index);

% Find indeces of hull who satisfy
[maybe_po, final_pos] = deal(cell(1, index));
for i = 1:index
    if ~isempty(hulls{1, i})
        maybe_po{i} = find(lbound{i} - ubound{i} <= 0);
    end
end

if minval ~= 0
    for i = 1:index
        if ~isempty(maybe_po{i})
            po = (minval - hulls{3, i}(maybe_po{i}))./abs(minval) +...
                hulls{2, i}(maybe_po{i}).*ubound{i}...
                (maybe_po{i})./abs(minval) >= SS.ep;
            final_pos{i} = hulls{1, i}(po);
        end
    end
else
    for i = 1:index
        if ~isempty(maybe_po{i})
            po = (hulls{3, i}(maybe_po{i})) -...
                (hulls{2, i}(maybe_po{i})).*ubound{i}(maybe_po{i}) <= 0;
            final_pos{i} = hulls{1, i}(po);
        end
    end
end
return

%--------------------------------------------------------------------------
% Function   :  calc_ubound
% Purpose    :  calculate the ubound used in determing potentially
%               optimal hrectangles
%--------------------------------------------------------------------------
function ub = calc_ubound(fc, hull, szes, index)
hull_lengths = szes;
ub           = cell(1, index);

for i = 1:index
    if ~isempty(hull{1, i})
        for j = 1:size(hull{1, i}, 2)
            tmp_rects = find(hull_lengths > hull{2, i}(j));
            if ~isempty(tmp_rects)
                tmp_f    = fc(tmp_rects);
                tmp_szes = szes(tmp_rects);
                tmp_ubs  = (tmp_f - hull{3, i}(j))./(tmp_szes -...
                            hull{2, i}(j));
                ub{i}(j) = min(tmp_ubs);
            else
                ub{i}(j) = 1.976e14;
            end
        end
    end
end
return

%--------------------------------------------------------------------------
% Function   :  calc_lbound
% Purpose    :  calculate the lbound used in determing potentially
%               optimal hrectangles
%--------------------------------------------------------------------------
function lb = calc_lbound(fc, hull, szes, index)
hull_lengths = szes;
lb           = cell(1, index);

for i = 1:index
    if ~isempty(hull{1, i})
        for j = 1:size(hull{1, i}, 2)
            tmp_rects = find(hull_lengths < hull{2, i}(j));
            if ~isempty(tmp_rects)
                tmp_f    = fc(tmp_rects);
                tmp_szes = szes(tmp_rects);
                tmp_ubs  = (hull{3, i}(j) - tmp_f)./(hull{2, i}(j) -...
                            tmp_szes);
                lb{i}(j) = max(tmp_ubs);
            else
                lb{i}(j) = -1.976e14;
            end
        end
    end
end
return

%--------------------------------------------------------------------------
% Function  : Evaluations
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, S, VAL] = Evaluations(Const, i, MSS, third, VAL, PH, OPTI)
%--------------------------------------------------------------------------
% Create and allocating empty sets
A.cO = MSS(i).C(:, PH);
A.LO = MSS(i).L(:, PH);
A.eO = MSS(i).E(PH);
A.O = MSS(i).F(PH);
A.Ox = MSS(i).X(PH);
DIVIS = find(A.LO == min(A.LO(:,1)));
DELTA = third(min(A.LO) + 1);
S = length(DIVIS);
[A.L, A.R, A.Lx, A.Rx, A.eL, A.eR] = deal(zeros(1, S));
A.lL = A.LO(:, 1)*ones(1, S);
A.cL = A.cO(:, 1)*ones(1, S); A.cR = A.cL;

for g = 1:S
    A.cL(DIVIS(g), g) = A.cL(DIVIS(g), g) - DELTA;
    om_point = abs(VAL.b - VAL.a).*A.cL(:, g) + VAL.a;
    if isfield(Const, 'constraint')
        A.Lx(g) = CallC(Const, om_point, OPTI);
        if A.Lx(g) == 0
            A.L(g) = feval(Const.f, om_point);
        else
            A.Lx(g) = 1;
            A.L(g) = 10^9;
        end
    else
        A.L(g) = feval(Const.f, om_point);
        if ~isnan(A.L(g)) && ~isinf(A.L(g)) 
            A.Lx(g) = 0;
        else
            A.Lx(g) = 1;
            A.L(g) = 10^9;
        end
    end
    
    A.cR(DIVIS(g), g) = A.cR(DIVIS(g), g) + DELTA;
    om_point = abs(VAL.b - VAL.a).*A.cR(:, g) + VAL.a;
    if isfield(Const, 'constraint')
        A.Rx(g) = CallC(Const, om_point, OPTI);
        if A.Rx(g) == 0
            A.R(g) = feval(Const.f, om_point);
        else
            A.Rx(g) = 1;
            A.R(g) = 10^9;
        end
    else
        A.R(g) = feval(Const.f, om_point);
        if ~isnan(A.R(g)) && ~isinf(A.R(g)) 
            A.Rx(g) = 0;
        else
            A.Rx(g) = 1;
            A.R(g) = 10^9;
        end
    end 
end
[~, A.lsx] = sort([min(A.L, A.R)' DIVIS], 1);
for g = 1:S
    A.lL(DIVIS(A.lsx(1:g, 1)), g) = A.lL(DIVIS(A.lsx(1:g, 1)), g) + 1;
    A.eL(g) = VAL.fcount + 1;
    A.eR(g) = VAL.fcount + 2;
    VAL.fcount = A.eR(g);
end
A.LO = A.lL(:, S);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  CallConstraints
% Purpose    :  Evaluate Constraints at pointed specified
%--------------------------------------------------------------------------
function ff_value = CallC(Problem, point, OPTI)
%--------------------------------------------------------------------------
ret_value = 0;
if isfield(Problem, 'constraint')
    [con_g, con_h] = feval(Problem.constraint, point);
    % inequality
    for i = 1:length(con_g)
        if con_g(i) > OPTI.ept
            ret_value = ret_value + con_g(i);
        end
    end
    % equality
    for i = 1:length(con_h)
        if abs(con_h(i)) > OPTI.ept
            ret_value = ret_value + abs(con_h(i));
        end
    end
end
if ret_value == 0
    ff_value = 0;
else
    ff_value = 1;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_min
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function [MV, minval, xatmin] = Find_min(CE, MSS)
%--------------------------------------------------------------------------
MV = nan(2, size(MSS, 2));

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        [MV(1, i), MV(2, i)] = min(MSS(i).F(1:CE(i)));
    end
end

[minval, Least] = min(MV(1, :));
xatmin = MSS(Least).C(:, MV(2, Least));
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  replaceinf                                         
% Purpose    :  Assign R. Carter value to given point              
%--------------------------------------------------------------------------
function MSS = replaceinf(MSS, CE, VAL, pert)
% Initialize fcn_values to original values
feas_pointsf = cell(1, size(MSS, 2));
feas_points  = cell(1, size(MSS, 2));
for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        XX = find(MSS(i).X(1:CE(i)) == 0);
        feas_points{i}  = MSS(i).C(:, XX);
        feas_pointsf{i} = MSS(i).F(XX);
    end
end

feas_points   = cell2mat(feas_points);
feas_pointsf = cell2mat(feas_pointsf);

if ~isempty(feas_pointsf)
    maxfc = max(feas_pointsf);
else
    maxfc = 10^9;
end

if ~isempty(feas_points)
    for i = 1:size(MSS, 2)
        if CE(i) ~= 0
            XX = find(MSS(i).X(1:CE(i)) == 1);
            for h = 1:length(XX)
                index = XX(h);
                found_pp = feas_points;
                found_ff = feas_pointsf;
                for j = 1:VAL.n
                    neighbors = find(abs(found_pp(j, :) - MSS(i).C(j,...
                    index)) <= 3^(-MSS(i).L(j, index)));
                    if ~isempty(neighbors)
                        found_pp  = found_pp(:, neighbors);
                        found_ff = found_ff(neighbors);
                    else
                        found_ff = [];
                        break;
                    end
                end
                
                % Assign Carter value to the point
                if ~isempty(found_ff)
                    % assign to index the min. value found + a little bit more
                    fstar = min(found_ff);
                    if fstar ~= 0
                        MSS(i).F(index) = fstar + pert*abs(fstar);
                    else
                        MSS(i).F(index) = fstar + pert*1;
                    end
                else
                    MSS(i).F(index)   = maxfc + 1;
                    maxfc             = maxfc + 1;
                end
            end
        end
    end
else
    for i = 1:size(MSS, 2)
        if CE(i) ~= 0
            XX = find(MSS(i).X(1:CE(i)) == 1);
            for h = 1:length(XX)
                MSS(i).F(XX(h))   = maxfc + 1;
                maxfc             = maxfc + 1;
            end
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% GETOPTS Returns options values in an options structure
%--------------------------------------------------------------------------
function varargout = getopts(options, varargin)
K = fix(nargin/2);
if nargin/2 == K
    error('fields and default values must come in pairs')
end
if isa(options,'struct')
    optstruct = 1;
else
    optstruct = 0;
end
varargout = cell(K, 1);
k = 0; ii = 1;
for i = 1:K
    if optstruct && isfield(options, varargin{ii})
        assignin('caller', varargin{ii}, options.(varargin{ii}));
        k = k + 1;
    else
        assignin('caller', varargin{ii}, varargin{ii + 1});
    end
    ii = ii + 2;
end
if optstruct && k ~= size(fieldnames(options), 1)
    warning('options variable contains improper fields')
end
return
%--------------------------------------------------------------------------
% END of BLOCK
%--------------------------------------------------------------------------