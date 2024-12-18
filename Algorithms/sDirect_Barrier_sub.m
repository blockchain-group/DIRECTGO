function [minima, xatmin, history] = sDirect_Barrier_sub...
    (Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : sDirect_Barrier_sub
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for problems with various
%              constraints constraints
%--------------------------------------------------------------------------
% [minima, xatmin, history] = sDirect_Barrier_sub(Problem, opts, bounds)
%       s         - static memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       Barrier   - Barrier approach for problems with constraints
%       sub       - Sub-dividing Step
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
%                 opts.showits   = 1 print iteration status
%                 opts.ep        = global/local weight parameter
%                 opts.ept       = tollerance for constrains
%                 opts.tol       = tolerance for termination if testflag=1               
%                 opts.sub       = set of iterations, in which subdividing
%                 step will be performed [it_1, it_2,...,it_k]
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
%
% J. Na, Y. Lim, C. Han. "A modi?ed DIRECT algorithm for hidden constraints
% in an LNG process optimization". Energy (2017) 488�500
% doi:10.1016/j.energy.2017.03.047.
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[third, VAL, MSS] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS] = Initialization(VAL, OPTI,...
    Problem, MSS);

while VAL.perror > OPTI.TOL                                 % Main loop
%--------------------------------------------------------------------------
    % Selection of potential optimal hyper-rectangles step
    POH = Find_po(MSS.FF(1:VAL.I), MSS.LL(:,1:VAL.I), Fmin, OPTI.ep,...
        MSS.DD(1:VAL.I));
    if ~isempty(OPTI.sub)
        if VAL.itctr == OPTI.sub(1)
            II = find(MSS.XX(1:VAL.I) == 1);
            POH = union(POH, II);
            OPTI.sub(1) = [];
        end
    end
    
    % Subdivide potential optimalhyper-rectangles
    [MSS,VAL] = Subdivision(VAL, Problem, third, MSS, POH, OPTI);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS);
    
%--------------------------------------------------------------------------
end                                                         % End of while

% Return value
minima  = Fmin;
if OPTI.G_nargout == 2
    xatmin  = Xmin;
elseif OPTI.G_nargout == 3
    xatmin  = Xmin;
    history = VAL.history(1:(VAL.itctr), 1:4);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Options
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
    1e-8, 'globalmin', 0, 'globalxmin', 0, 'ep', 1e-4, 'sub',...
    [5, 25, 625, 3125, 15625, 78125, 390625, 1953125]);

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
        OPTI.globalXMIN = getInfo.xmin(VAL.n);
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
OPTI.sub         = sub;    % numbers of iterations in which sub
                           % dividing step will be performed
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : BEGIN
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [third, VAL, MSS] = Alocate(OPTI, VAL)
%--------------------------------------------------------------------------
tic                                     % Mesure time
VAL.time    = 0;                        % initial time
VAL.itctr   = 1;                        % initial iteration
VAL.perror  = 10;                       % initial perror
VAL.barrier = 10^9;                     % initial barrier value

% alociate MAIN sets
z   = round(OPTI.MAXevals);
MSS = struct('FF', zeros(1, z), 'XX', zeros(1, z), 'II', zeros(1, z),...
    'DD', -ones(1, z), 'LL', zeros(VAL.n, z), 'CC', zeros(VAL.n, z));

third       = zeros(1, 10*VAL.n*OPTI.MAXdeep);   % delta values
third(1)    = 1/3;                               % first delta
for i = 2:10*VAL.n*OPTI.MAXdeep                  % all delta
    third(i)  = (1/3)*third(i - 1);
end
if OPTI.G_nargout == 3
    VAL.history = zeros(OPTI.MAXits, 4);  % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [OPTI, VAL, Xmin, Fmin, MSS] = Initialization(VAL, OPTI,...
    Problem, MSS)
%--------------------------------------------------------------------------
VAL.I          = 1;                     % evaluation counter
MSS.DD(1)      = 1;
MSS.II(1)      = 1;
MSS.CC(:, 1)   = ones(VAL.n, 1)/2;
Point = abs(VAL.b - VAL.a).*(MSS.CC(:, 1)) + VAL.a;
if isfield(Problem, 'constraint')
    MSS.XX(1) = CallC(Problem, Point, OPTI);
    if MSS.XX(1) == 0
        MSS(1).FF(1) = deal(feval(Problem.f, Point));
    else
        MSS(1).FF(1) = deal(10^9);
    end
else
    ff = feval(Problem.f, Point);
    if ~isnan(ff) && ~isinf(ff) 
        MSS(1).FF(1) = deal(ff);
        MSS.XX(1) = 0;
    else
        MSS.XX(1) = 1;
        MSS(1).FF(1) = deal(10^9);
    end
end

[VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS);

% Check stop condition if global minima is known
if OPTI.TESTflag  == 1
    if OPTI.globalMIN ~= 0
        VAL.perror = 100*(Fmin - OPTI.globalMIN)/abs(OPTI.globalMIN);
    else
        VAL.perror = 100*Fmin;
    end
else
    VAL.perror   = 2;
end
if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.I;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
function [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS)
%--------------------------------------------------------------------------
[Fmin, fminindex] = min(MSS.FF(1:VAL.I));
Xmin = (abs(VAL.b - VAL.a)).*MSS.CC(:, fminindex) + VAL.a;
%--------------------------------------------------------------------------
VAL.time = toc;

if OPTI.showITS == 1               % Show iteration stats
    fprintf(...
        'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, Fmin, VAL.time, VAL.I);
end
if OPTI.TESTflag == 1               % Check for stop condition
    if OPTI.globalMIN ~= 0            % Calculate error if globalmin known
        VAL.perror = 100*(Fmin - OPTI.globalMIN)/abs(OPTI.globalMIN);
    else
        VAL.perror = 100*Fmin;
    end
    if VAL.perror < OPTI.TOL
        fprintf('Minima was found with Tolerance: %4i', OPTI.TOL);
        VAL.perror = -10;
    end
else
    VAL.perror = 10;
end

if VAL.itctr >= OPTI.MAXits         % Have we exceeded the maxits?
    disp('Exceeded max iterations. Increase maxits');
    VAL.perror = -10;
end

if VAL.I > OPTI.MAXevals       % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if max(max(MSS.LL)) > OPTI.MAXdeep  % Have we exceeded max deep?
    disp('Exceeded Max depth. Increse maxdeep');
    VAL.perror = -10;
end

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.I;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end

% Update iteration number
if VAL.perror > OPTI.TOL
    VAL.itctr = VAL.itctr + 1;
end
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
% Function  : Subdivision
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [MSS, VAL] = Subdivision(VAL, O, third, MSS, POHas, OPTI)
%--------------------------------------------------------------------------
POH           = size(POHas, 2);
[c_lt, c_rt, f_l, f_r, col, cor, order, s_tp] = deal(cell(1, POH));
li            = (MSS.LL(:, POHas));
oldc          = (MSS.CC(:, POHas));
[l_tp, ls]    = deal(cell(1, POH));
[Delta, LL]   = deal(zeros(POH, 1));
FCX           = zeros(POH + 1, 1);
for i = 1:POH
    ls{i}       = find(li(:, i) == min(li(:, i)));
    LL(i)       = length(ls{i});
    FCX(i)      = VAL.I + sum(LL(1:(i - 1)))*2;
    l_tp{i}     = li(:, i)*ones(1, LL(i));
    Delta(i)    = third(min(li(: ,i)) + 1);
    [f_l{i}, f_r{i}, col{i}, cor{i}] = deal(zeros(1, LL(i)));
end
FCX(POH + 1)  = VAL.I + sum(LL)*2;
VAL.I         = FCX(end);

% Calculate new points and evaluate at the objective function
for i = 1:POH
    dd                = diag(repelem(Delta(i), LL(i)));
    c_lt{i}           = oldc(:, i)*ones(1, LL(i));
    c_rt{i}           = c_lt{i};
    c_lt{i}(ls{i}, :) = c_lt{i}(ls{i}, :) - dd;
    c_rt{i}(ls{i}, :) = c_rt{i}(ls{i}, :) + dd;
    point_l           = abs(VAL.b - VAL.a).*c_lt{i} + VAL.a;
    point_r           = abs(VAL.b - VAL.a).*c_rt{i} + VAL.a;
    
    for j = 1:LL(i)
        if isfield(O, 'constraint')
            col{i}(j) = CallC(O, point_l(:, j), OPTI);
            if col{i}(j) == 0
                f_l{i}(j) = feval(O.f, point_l(:, j));
            else
                col{i}(j) = 1;
                f_l{i}(j) = VAL.barrier;
                VAL.barrier = VAL.barrier + 1;
            end
        else
            f_l{i}(j) = feval(O.f, point_l(:, j));
            if ~isnan(f_l{i}(j)) && ~isinf(f_l{i}(j))
                col{i}(j) = 0;
            else
                col{i}(j) = 1;
                f_l{i}(j) = VAL.barrier;
                VAL.barrier = VAL.barrier + 1;
            end
        end
        if isfield(O, 'constraint')
            cor{i}(j) = CallC(O, point_r(:, j), OPTI);
            if cor{i}(j) == 0
                f_r{i}(j) = feval(O.f, point_r(:, j));
            else
                cor{i}(j) = 1;
                f_r{i}(j) = VAL.barrier;
                VAL.barrier = VAL.barrier + 1;
            end
        else
            f_r{i}(j) = feval(O.f, point_r(:, j));
            if ~isnan(f_r{i}(j)) && ~isinf(f_r{i}(j))
                cor{i}(j) = 0;
            else
                cor{i}(j) = 1;
                f_r{i}(j) = VAL.barrier;
                VAL.barrier = VAL.barrier + 1;
            end
        end
    end
    
    [~, order{i}] = sort([min(((f_l{i})), ((f_r{i})))' ls{i}], 1);
    for j = 1:LL(i)
        l_tp{i}(ls{i}(order{i}(1:j, 1)), j) = l_tp{i}(ls{i}(order{i}...
            (1:j, 1)), j) + 1;
        s_tp{i}(j) = 1/2*norm((1/3*(ones(VAL.n, 1))).^(l_tp{i}(:, j)));
    end
end

% Store information
for i=1:POH
    f_tp = f_r{i}(:, [1; 1]*transpose(order{i}(:, 1)));
    f_tp(:, 1:2:end) = f_l{i}(:, order{i}(:, 1));
    MSS.FF(FCX(i) + 1:FCX(i + 1)) = f_tp;                  % f(x)
    x_tp = cor{i}(:,[1; 1]*transpose(order{i}(:, 1)));
    x_tp(:, 1:2:end) = col{i}(:, order{i}(:, 1));
    MSS.XX(FCX(i) + 1:FCX(i + 1)) = x_tp;                  % g(x)
    c_tp = c_rt{i}(:, [1; 1]*transpose(order{i}(:, 1)));
    c_tp(:, 1:2:end) = c_lt{i}(:, order{i}(:, 1));
    MSS.CC(:, FCX(i) + 1:FCX(i + 1)) = c_tp;               % x
    MSS.LL(:, FCX(i) + 1:FCX(i + 1)) = repelem(l_tp{i}, 1, 2);
    MSS.LL(:, POHas(i)) = l_tp{i}(:, LL(i));               % Side lengths
    MSS.DD(FCX(i) + 1:FCX(i + 1)) = repelem(s_tp{i}, 1, 2);
    MSS.DD(POHas(i)) = s_tp{i}(LL(i));                     % Diameterds
    MSS.II(FCX(i) + 1:FCX(i + 1)) = FCX(i) + 1:FCX(i + 1); % indexes
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_po                                            
% Purpose    :  Return list of PO hyperrectangles                 
%--------------------------------------------------------------------------
function final_pos = Find_po(fc, lengths, minval, ep, szes)
% Find all rects on hub
diff_szes = sum(lengths, 1);
tmp_max = max(diff_szes);
j = 1;
sum_lengths = sum(lengths, 1);
hull = zeros(1, size(fc, 2));
for i = 1:tmp_max + 1
    tmp_idx = find(sum_lengths == i - 1);
    [tmp_n, hullidx] = min(fc(tmp_idx));
    if ~isempty(hullidx)
        hull(j) = tmp_idx(hullidx);
        j = j + 1;
        % Check for ties
        ties = find(abs(fc(tmp_idx) - tmp_n) <= 1e-13);
        if length(ties) > 1
            mod_ties = find(tmp_idx(ties) ~= hull(j - 1));
            hull(j: j + size(mod_ties, 2) - 1) = tmp_idx(ties(mod_ties));
            j = j + size(mod_ties,2);
        end
    end
end
hull = hull(1:j - 1);

% Compute lb and ub for rects on hub
lbound = calc_lbound(lengths, fc, hull, szes);
ubound = calc_ubound(lengths, fc, hull, szes);

% Find indeces of hull who satisfy
maybe_po = find(lbound-ubound <= 0);

if minval ~= 0
    po = find((minval - fc(hull(maybe_po)))./abs(minval) +...
        szes(hull(maybe_po)).*ubound(maybe_po)./abs(minval) >= ep);
else
    po = find(fc(hull(maybe_po)) -...
        szes(hull(maybe_po)).*ubound(maybe_po) <= 0);
end
final_pos      = hull(maybe_po(po));
return

%--------------------------------------------------------------------------
% Function   :  calc_ubound                                                                          
% Purpose    :  calculate the ubound used in determing potentially 
%               optimal hrectangles                               
%--------------------------------------------------------------------------
function ub = calc_ubound(lengths, fc, hull, szes)
hull_length  = length(hull);
hull_lengths = lengths(:, hull);
ub           = zeros(1, hull_length);
for i =1:hull_length
    tmp_rects = find(sum(hull_lengths, 1) < sum(lengths(:, hull(i))));
    if ~isempty(tmp_rects)
        tmp_f     = fc(hull(tmp_rects));
        tmp_szes  = szes(hull(tmp_rects));
        tmp_ubs   = (tmp_f - fc(hull(i)))./(tmp_szes - szes(hull(i)));
        ub(i)     = min(tmp_ubs);
    else
        ub(i)     = 1.976e14;
    end
end
return

%--------------------------------------------------------------------------
% Function   :  calc_lbound                                                                            
% Purpose    :  calculate the lbound used in determing potentially 
%               optimal hrectangles                               
%--------------------------------------------------------------------------
function lb = calc_lbound(lengths, fc, hull, szes)
hull_length  = length(hull);
hull_lengths = lengths(:, hull);
lb           = zeros(1, hull_length);
for i = 1:hull_length
    tmp_rects = find(sum(hull_lengths,1) > sum(lengths(:, hull(i))));
    if ~isempty(tmp_rects)
        tmp_f     = fc(hull(tmp_rects));
        tmp_szes  = szes(hull(tmp_rects));
        tmp_lbs   = (fc(hull(i)) - tmp_f)./(szes(hull(i)) - tmp_szes);
        lb(i)     = max(tmp_lbs);
    else
        lb(i)     = -1.976e14;
    end
end
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