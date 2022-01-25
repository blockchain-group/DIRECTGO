function [minima, xatmin, history] = sDirect_GLce_min(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : sDirect_GLce_min
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for problems with various
%              constraints constraints
%--------------------------------------------------------------------------
% [minima, xatmin, history] = sDirect_GLce_min(Problem, opts, bounds)
%       s         - static memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       G         - Enhancing the global search
%       L         - Enhancing the local search
%       c         - Constraint handling
%       e         - Additional tolerance for constraint handling
%       min       - Matlab optimization solver fmincon
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
%                 opts.showits   = 1 print iteration status
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
%
% Constraint handling strategie taken from:
%--------------------------------------------------------------------------
% Stripinis, L., Paulavicius, R., Zilinskas, J.: "Penalty functions and
% two-step selection procedure based DIRECT-type algorithm for constrained
% global optimization".  Structural and Multidisciplinary Optimization,
% (2019). ISSN 1615-1488, DOI: 10.1007/s00158-018-2181-2
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[third, VAL, MSS] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS, XFmin] = Initialization(VAL, OPTI, Problem,...
    MSS, third);

while VAL.perror > OPTI.TOL                                 % Main loop
%--------------------------------------------------------------------------
    % Selection of potential optimal hyper-rectangles step
    POH = Selection(VAL, MSS, Xmin);
    
    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL] = Subdivision(VAL, Problem, third, MSS, POH, OPTI, 2);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin, MSS, XFmin] = Arewedone(OPTI, VAL, MSS, Problem);
    
%--------------------------------------------------------------------------
end                                                         % End of while

% Return value
minima      = Fmin;
if OPTI.G_nargout == 2
    xatmin    = XFmin;
elseif OPTI.G_nargout == 3
    xatmin    = XFmin;
    history   = VAL.history(1:(VAL.itctr), 1:4);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% AUXILIARY FUNCTION BLOCK
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Function  : OPTI
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
    0, 'globalmin', 0, 'globalxmin', 0);

if isempty(bounds)

% Return the problem information.
getInfo = feval(Problem.f);

if getInfo.ng ~= 0 || getInfo.nh ~= 0
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
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [third, VAL, MSS] = Alocate(OPTI, VAL)
%--------------------------------------------------------------------------
tic                                     % Mesure time
VAL.time    = 0;                        % initial time
VAL.itctr   = 1;                        % initial iteration
VAL.perror  = 10;                       % initial perror
VAL.epsil   = 1;                        % initial espil {0.0001<epsil<1}
VAL.ep      = 0.0001;
VAL.CARD    = 10*(VAL.n^3);             % infeas HR
VAL.CARDi   = 1000*(VAL.n^3);
VAL.tet     = 10^(-6);
VAL.STAG    = 0;                        % iteration stagnates?

% alociate MAIN sets
z   = round(OPTI.MAXevals);
MSS = struct('FF', zeros(1, z), 'df', zeros(1, z), 'XX', zeros(1, z),...
    'II', zeros(1, z), 'DD', -ones(1, z), 'LL', zeros(VAL.n, z),...
    'CC', zeros(VAL.n, z));

third       = zeros(1, 10*VAL.n*OPTI.MAXdeep);   % delta values
third(1)    = 1/3;                               % first delta
for i = 2:10*VAL.n*OPTI.MAXdeep                  % all delta
    third(i)  = (1/3)*third(i - 1);
end
if OPTI.G_nargout == 3
    VAL.history = zeros(OPTI.MAXits, 4); % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Selection
% Purpose   : Selection of potential optimal hyper-rectangles
%--------------------------------------------------------------------------
function POH = Selection(VAL, MSS, Xmin)
%--------------------------------------------------------------------------
% Calculate Euclidean Distatnces
Euclid_dist = sum((Xmin - MSS.CC(:, 1:VAL.I)).^2, 1).^0.5;

% Identify potential optimal hyper-rectangles
S = Find_poh(MSS.FF(1:VAL.I) + MSS.df(1:VAL.I), MSS.DD(1:VAL.I),...
    MSS.II(1:VAL.I));

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
function [OPTI, VAL, Xmin, Fmin, MSS, XFmin] = Initialization(VAL, OPTI,...
    Problem, MSS, third)
%--------------------------------------------------------------------------
VAL.I          = 1;                       % evaluation counter
MSS.DD(1)      = 1;                       % initial diameter
MSS.II(1)      = 1;                       % initial index
MSS.df(1)      = 0;
MSS.CC(:, 1)   = ones(VAL.n, 1)/2;        % initial midpoint
Point          = abs(VAL.b - VAL.a).*(MSS.CC(:, 1)) + VAL.a;  
[MSS(1).FF(1), VAL.min_last, VAL.fminloc,...
    VAL.min_gen] = deal(feval(Problem.f, Point));
MSS.XX(1)      = CallC(Problem, Point, OPTI);
Xmin           = MSS.CC(:, 1);            % initial point
VAL.stagnate = 0;
VAL.nFunc = 0;
VAL.Xminval = Point;
VAL.local = 0;
if MSS.XX(1) ~= 0
    id = 0;
    fprintf('Phase II: searching feasible point:'); disp(' ');
    while id == 0
%--------------------------------------------------------------------------
        % Calculate Euclidean Distatnces
        Euclid_dist = sum((Xmin - MSS.CC(:, 1:VAL.I)).^2, 1).^0.5;
        
        % Identify potential optimal hyper-rectangles using normalized values
        S = Find_poh(MSS.XX(1:VAL.I), MSS.DD(1:VAL.I), MSS.II(1:VAL.I));
        
        % Identify potential optimal hyper-rectangles using Euclidean distances
        D = Find_poh(Euclid_dist(1:VAL.I), MSS.DD(1:VAL.I), MSS.II(1:VAL.I));
        
        % Find unique set of potential optimal hyper-rectangles
        POH = union(S, D);
%--------------------------------------------------------------------------
        [MSS, VAL] = Subdivision(VAL, Problem, third, MSS, POH, OPTI, 1);
        
        Fmin          = min(MSS.XX(1:VAL.I));
        fminindex     = find(MSS.XX(1:VAL.I) == Fmin);
        if size(fminindex, 2) ~= 1
            I           = min(MSS.DD(fminindex));
            TI          = MSS.DD(fminindex) == I;
            fminindex   = fminindex(TI);
            if size(fminindex, 2) ~= 1
                [~, I]    = max((MSS.II(fminindex)));
                fminindex = fminindex(I);
            end
        end
        Xmin          = MSS.CC(:, fminindex);
        
        if Fmin == 0
            VAL.local = 1;
            Fmin = min(MSS.FF(MSS.XX(1:VAL.I) == 0));          
            id = 1;
            fprintf('Violation of constraints: %15.10f  feasible point: %15.10f  fn evals: %8i\n',...
                0, Fmin, VAL.I);
        else
            fprintf('Violation of constraints: %15.10f    fn evals: %8i\n',...
                Fmin, VAL.I);
        end
    end
    fprintf('Phase I: Improve feasible solution:'); disp(' ');
end
VAL.Xmin = Xmin;
VAL.Xminval = abs(VAL.b - VAL.a).*Xmin + VAL.a;  
[VAL.min_last, VAL.fminloc, VAL.min_gen] = deal(Fmin);
% Check stop condition if global minima is known
[VAL, Fmin, Xmin, MSS, XFmin] = Arewedone(OPTI, VAL, MSS, Problem);

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.I + VAL.nFunc;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
function [VAL, Fmin, Xmin, MSS, XFmin] = Arewedone(OPTI, VAL, MSS, Problem)
%--------------------------------------------------------------------------
TMf           = find(MSS.XX(1:VAL.I) == 0);
Fmin          = min(MSS.FF(TMf));
index         = TMf(find(MSS.FF(TMf) == Fmin, 1, 'last'));
XFmin         = abs(VAL.b - VAL.a).*(MSS.CC(:, index)) + VAL.a;  
MSS.df        = zeros(1, VAL.I);
TMe = find(MSS.XX(1:VAL.I) > 0 & MSS.XX(1:VAL.I) < VAL.epsil);

if VAL.STAG == 10 && VAL.CARD <= VAL.CARDi
    VAL.STAG    = 0;
    VAL.epsil   = 1;
    VAL.CARD    = VAL.CARD * 10;
    VAL.tet     = VAL.tet/100;
elseif (size(TMe, 2) == 0)       && (VAL.epsil * 3 <= 10)
    VAL.epsil   = VAL.epsil * 3;
elseif (size(TMe, 2) > VAL.CARD) && (VAL.epsil / 3 >= VAL.ep)
    VAL.epsil   = VAL.epsil / 3;
elseif (size(TMe, 2) > VAL.CARD) && (VAL.epsil / 3 <= VAL.ep)
    VAL.epsil = 0;
end
TM = union(TMf, TMe);
TMc = setdiff(1:VAL.I, TM);

MSS.df(TMc) = abs(MSS.FF(TMc) - Fmin) + MSS.XX(TMc);
Least = find((MSS.FF(1:VAL.I) + MSS.df(1:VAL.I)) ==...
            min(MSS.FF(1:VAL.I) + MSS.df(1:VAL.I)), 1, 'last');
Xmin = MSS.CC(:, Least);

if VAL.epsil == 0
    LEMDA = abs(((abs(VAL.b - VAL.a).*Xmin + VAL.a))...
            - (abs(VAL.b - VAL.a).*VAL.Xmin + VAL.a));
    if sum((LEMDA).^2)^0.5 < VAL.tet
        VAL.STAG = VAL.STAG + 1;
    else
        VAL.STAG = 0;
    end
end
VAL.Xmin = Xmin;

if (abs(Fmin - VAL.min_last) < 0.01*abs(Fmin))
    VAL.stagnate = VAL.stagnate + 1;
else
    VAL.stagnate    = 0;
    VAL.min_last = Fmin;
end

[Fmin, VAL.min_gen] = deal(min(Fmin, VAL.min_gen));
if Fmin < VAL.fminloc
    VAL.Xminval = XFmin; 
end
if (VAL.stagnate > 0) && (abs(VAL.min_gen - VAL.fminloc) >...
        0.01*abs(VAL.min_gen)) || VAL.local == 1
    VAL.local = 0;
    
    options = optimoptions('fmincon', 'Display', 'none', 'ConstraintTolerance', OPTI.ept);
    [xminloc, VAL.fminloc, ~, output] = fmincon(Problem.f, XFmin, [],...
        [], [], [], VAL.a, VAL.b, Problem.constraint, options);
    
    VAL.nFunc = VAL.nFunc + output.funcCount;
    if (VAL.fminloc < VAL.min_gen)
        VAL.min_gen = VAL.fminloc;
        Fmin = VAL.fminloc;
        XFmin = xminloc; 
        VAL.Xminval = xminloc;
    end
end
%--------------------------------------------------------------------------
if OPTI.showITS == 1               % Show iteration stats
    VAL.time = toc;
    fprintf(...
        'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, Fmin, VAL.time, VAL.I + VAL.nFunc);
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
    VAL.history(VAL.itctr, 2) = VAL.I + VAL.nFunc;
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
function [MSS, VAL] = Subdivision(VAL, O, third, MSS, POHas, OPTI, P)
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
end
FCX(POH + 1)  = VAL.I + sum(LL)*2;
VAL.I         = FCX(end);

% Calculate new points and evaluate at the objective function
for i = 1:POH
    dd = diag(repelem(Delta(i), LL(i)));
    c_lt{i} = oldc(:, i)*ones(1, LL(i)); c_rt{i} = c_lt{i};
    c_lt{i}(ls{i}, :) = c_lt{i}(ls{i},:) - dd;
    c_rt{i}(ls{i}, :) = c_rt{i}(ls{i},:) + dd;
    point_l = abs(VAL.b - VAL.a).*c_lt{i} + VAL.a;
    point_r = abs(VAL.b - VAL.a).*c_rt{i} + VAL.a;
    col{i} = transpose(arrayfun(@(x) CallC(O, point_l(:,x),...
        OPTI), (1:LL(i)).'));
    cor{i} = transpose(arrayfun(@(x) CallC(O, point_r(:,x),...
        OPTI), (1:LL(i)).'));
    f_l{i} = transpose(arrayfun(@(x)...
        feval(O.f, point_l(:,x)), (1:LL(i)).'));
    f_r{i} = transpose(arrayfun(@(x)...
        feval(O.f, point_r(:,x)), (1:LL(i)).'));
    if P == 2
        [~, order{i}] = sort([min(f_l{i}, f_r{i})' ls{i}], 1);
    else
        [~, order{i}] = sort([min(col{i}, cor{i})' ls{i}], 1);
    end
    for j = 1:LL(i)
        l_tp{i}(ls{i}(order{i}(1:j,1)), j) = l_tp{i}(ls{i}(order{i}...
            (1:j, 1)), j) + 1;
        s_tp{i}(j) = 1/2*norm((1/3*(ones(VAL.n, 1))).^(l_tp{i}(:,j)));
    end
end

% Store information
for i = 1:POH
    f_tp = f_r{i}(:,[1; 1]*transpose(order{i}(:, 1)));
    f_tp(:, 1:2:end) = f_l{i}(:, order{i}(:, 1));
    MSS.FF(FCX(i) + 1:FCX(i + 1)) = f_tp;                   % f(x)
    x_tp = cor{i}(:,[1; 1]*transpose(order{i}(:, 1)));
    x_tp(:, 1:2:end) = col{i}(:, order{i}(:, 1));
    MSS.XX(FCX(i) + 1:FCX(i + 1)) = x_tp;                   % g(x)
    c_tp = c_rt{i}(:,[1; 1]*transpose(order{i}(:, 1)));
    c_tp(:, 1:2:end) = c_lt{i}(:,order{i}(:, 1));
    MSS.CC(:,FCX(i) + 1:FCX(i + 1)) = c_tp;                 % x
    MSS.LL(:,FCX(i) + 1:FCX(i + 1)) = repelem(l_tp{i}, 1, 2);
    MSS.LL(:, POHas(i)) = l_tp{i}(:,LL(i));                 % Side lengths
    MSS.DD(FCX(i) + 1:FCX(i + 1)) = repelem(s_tp{i}, 1, 2);
    MSS.DD(POHas(i)) = s_tp{i}(LL(i));                      % Diameterds
    MSS.II(FCX(i) + 1:FCX(i + 1)) = FCX(i) + 1:FCX(i + 1);  % indexes
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  CallConstraints
% Purpose    :  Evaluate Constraints at pointed specified
%--------------------------------------------------------------------------
function ret_value = CallC(Problem, point, OPTI)
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
        if abs(con_h(i)) > 10^(-8)
            ret_value = ret_value + abs(con_h(i));
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function:   Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function boxes = Find_poh(fc, szes, indexs)
%--------------------------------------------------------------------------
C     = unique(szes);
SMS   = zeros(3, size(C, 2));
m_set = [fc; szes; indexs];
for i = 1: size(C, 2)                             % reduce m_set
    MB          = m_set(:, m_set(2, :) == C(i));
    [MV, MI]    = min(MB(1, :));
    TM          = find(MB(1, :) == MV);
    if size(TM, 2) ~= 1
        [~, MI]   = max(MB(3, TM));
        MI        = TM(MI);
    end
    SMS(:, i)   = MB(:, MI);
end
SMS           = flip(SMS,2);
% Create sets and define stop condition
fc_min    = SMS(1, :);
s_i       = 0;
setas     = zeros(1, size(fc_min, 2));
index     = size(fc_min, 2);
% Find index set of potential optimal hyper-rectangles
while index ~= 0
    [m_m, index]  = min(fc_min(1:index));
    if ~isnan(m_m)
        s_i        = s_i + 1;
        setas(s_i) = index;
    end
    index    = index - 1;
end
setas     = setas(1:s_i);
boxes     = SMS(3, setas);
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