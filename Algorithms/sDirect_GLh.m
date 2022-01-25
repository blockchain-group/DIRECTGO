function [minima, xatmin, history] = sDirect_GLh(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : sDirect_GLh
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for problems with various
%              constraints constraints
%--------------------------------------------------------------------------
% [minima, xatmin, history] = sDirect_GLh(Problem, opts, bounds)
%       s         - static memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       G         - Enhancing the global search
%       L         - Enhancing the local search
%       h         - Hidden constraint handling
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
%                 opts.ept       = tollerance for constrains
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
% Selection of potential optimal hyper-rectangles taken from:
%--------------------------------------------------------------------------
% Stripinis, L., Paulavicius, R., Zilinskas, J.: Improved scheme for
% selection of potentially optimal hyperrectangles in DIRECT. Optimization
% Letters (2018). ISSN 1862-4472, 12 (7), 1699-1712,
% DOI: 10.1007/s11590-017-1228-4
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[third, VAL, MSS] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, Fmax, MSS] = Initialization(VAL, OPTI,...
    Problem, MSS, third);

while VAL.perror > OPTI.TOL                                 % Main loop
%--------------------------------------------------------------------------
    % Selection of potential optimal hyper-rectangles step
    [POH, MSS] = Selection(VAL, MSS, Fmin, Fmax, Xmin);
    
    % Subdivide potential optimalhyper-rectangles
    [MSS,VAL] = Subdivision(VAL, Problem, third, MSS, POH, Xmin, Fmin,...
        2, OPTI);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Fmax, Xmin] = Arewedone(OPTI, VAL, MSS);
    
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
    1e-8, 'globalmin', 0, 'globalxmin', 0);

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
VAL.max     = sum((zeros(VAL.n, 1) - ones(VAL.n, 1)).^2, 1).^0.5;
VAL.perror  = 10;                       % initial perror

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
% Function  : Selection
% Purpose   : Selection of potential optimal hyper-rectangles
%--------------------------------------------------------------------------
function [POH, MSS] = Selection(VAL, MSS, Fmin, Fmax, Xmin)
%--------------------------------------------------------------------------
% Calculate Euclidean Distatnces
Euclid_dist = sum((Xmin - MSS.CC(:, 1:VAL.I)).^2, 1).^0.5;

% Hidden constraints replace
II = MSS.XX(1:VAL.I) == 1;
MSS.FF(II) = Fmin + Euclid_dist(II);

% Normalize Objective function values
TT = (MSS.FF(1:VAL.I) - Fmin + 10^(-16))/(Fmax - Fmin + 10^(-16));
TT(II) = Euclid_dist(II)/VAL.max;
TT = TT.*MSS.DD(1:VAL.I);
% Identify potential optimal hyper-rectangles using normalized values
S = Find_poh(TT(1:VAL.I), MSS.DD(1:VAL.I), MSS.II(1:VAL.I));

% Identify potential optimal hyper-rectangles using Euclidean distances
D = Find_poh(Euclid_dist(1:VAL.I), MSS.DD(1:VAL.I), MSS.II(1:VAL.I));

% Find unique set of potential optimal hyper-rectangles
POH = union(S, D);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [OPTI, VAL, Xmin, Fmin, Fmax, MSS] = Initialization(VAL, OPTI,...
    Problem, MSS, third)
%--------------------------------------------------------------------------
VAL.I = 1;                     % evaluation counter
MSS.DD(1) = 1;
MSS.II(1) = 1;
MSS.CC(:, 1) = ones(VAL.n, 1)/2;
Xmin = MSS.CC(:, 1);
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

if MSS.XX(1) == 1                         % initial midpoint feasible?
%--------------------------------------------------------------------------
    Fmin         = MSS.FF(1);             % initial minima
    LH_index     = 1;
    fprintf('Phase II: searching feasible point');
    while LH_index ~= 0
        MI       = max(MSS.II(MSS.DD(1:VAL.I) == max(MSS.DD(1:VAL.I))));
        [MSS, VAL] = Subdivision(VAL, Problem, third, MSS, MI,...
            Xmin, Fmin, 1, OPTI);
        LH_index   = min(MSS.XX(1:VAL.I));
        fprintf('fn evals: %8i\n', VAL.I);
    end
%--------------------------------------------------------------------------
end
[VAL, Fmin, Fmax, Xmin] = Arewedone(OPTI, VAL, MSS);

% Check stop condition if global minima is known
if OPTI.TESTflag  == 1
    if OPTI.globalMIN ~= 0
        VAL.perror = (Fmin - OPTI.globalMIN)/abs(OPTI.globalMIN);
    else
        VAL.perror = Fmin;
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
function [VAL, Fmin, Fmax, Xmin] = Arewedone(OPTI, VAL, MSS)
%--------------------------------------------------------------------------
index = find(MSS.XX(1:VAL.I) == 0);
Fmin = min(MSS.FF(index));         
Fmax = max(MSS.FF(1:VAL.I));
m_idx = MSS.FF(index) == Fmin;
m_idx = index(m_idx);
fminindex = m_idx(find(MSS.DD(m_idx) == max(MSS.DD(m_idx)), 1, 'last'));
Xmin = MSS.CC(:, fminindex);
%--------------------------------------------------------------------------
if OPTI.showITS == 1               % Show iteration stats
    VAL.time = toc;
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
% Function  : Subdivision
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [MSS, VAL] = Subdivision(VAL, O, third, MSS, POHas, Xmin,...
    Fmin, ph, OPTI)
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
                f_l{i}(j) = 10^9;
            end
            if col{i}(j) == 1 && ph == 2
                f_l{i}(j) = Fmin + sum((Xmin - point_l(:, j)).^2, 1).^0.5;
            end
        else
            f_l{i}(j) = feval(O.f, point_l(:, j));
            if ~isnan(f_l{i}(j)) && ~isinf(f_l{i}(j))
                col{i}(j) = 0;
            else
                col{i}(j) = 1;
                f_l{i}(j) = 10^9;
            end
            if col{i}(j) == 1 && ph == 2
                f_l{i}(j) = Fmin + sum((Xmin - point_l(:, j)).^2, 1).^0.5;
            end
        end
        if isfield(O, 'constraint')
            cor{i}(j) = CallC(O, point_r(:, j), OPTI);
            if cor{i}(j) == 0
                f_r{i}(j) = feval(O.f, point_r(:, j));
            else
                cor{i}(j) = 1;
                f_r{i}(j) = 10^9;
            end
            if cor{i}(j) == 1 && ph == 2
                f_r{i}(j) = Fmin + sum((Xmin - point_r(:, j)).^2, 1).^0.5;
            end
        else
            f_r{i}(j) = feval(O.f, point_r(:, j));
            if ~isnan(f_r{i}(j)) && ~isinf(f_r{i}(j))
                cor{i}(j) = 0;
            else
                cor{i}(j) = 1;
                f_r{i}(j) = 10^9;
            end
            if cor{i}(j) == 1 && ph == 2
                f_r{i}(j) = Fmin + sum((Xmin - point_r(:, j)).^2, 1).^0.5;
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
% Function:   Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function boxes = Find_poh(fc, szes, indexs)
%--------------------------------------------------------------------------
C = unique(szes);
SMS = zeros(3, size(C, 2));
m_set = [fc; szes; indexs];
for i = 1:size(C, 2)                             % reduce m_set
    MB = m_set(:, m_set(2, :) == C(i));
    TM = find(MB(1, :) == min(MB(1, :)));
    MI = TM(MB(3, TM) == max(MB(3, TM)));
    SMS(:, i) = MB(:, MI);
end
SMS = flip(SMS,2);
% Create sets and define stop condition
fc_min = SMS(1, :);
s_i = 0;
setas = zeros(1, size(fc_min, 2));
index = size(fc_min, 2);
% Find index set of potential optimal hyper-rectangles
while index ~= 0
    [m_m, index]  = min(fc_min(1:index));
    if ~isnan(m_m)
        s_i = s_i + 1;
        setas(s_i) = index;
    end
    index = index - 1;
end
setas = setas(1:s_i);
boxes = SMS(3, setas);
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