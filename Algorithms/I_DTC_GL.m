function [minima, xatmin, history] = I_DTC_GL(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : I_DTC_GL
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 03/17/2020
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = I_DTC_GL(Problem, opts, bounds)
%       I_DTC - Hyper-rectangular partitioning based on 1-Dimensional 
%               Trisection and objective function evaluations at Center 
%               points
%       G     - Enhancing the global search
%       L     - Enhancing the local search
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f       = Objective function handle
%
%       opts    - MATLAB structure which contains options.
%                 opts.maxevals  = max. number of function evals
%                 opts.maxits    = max. number of iterations
%                 opts.maxdeep   = max. number of rect. divisions
%                 opts.testflag  = 1 if globalmin known, 0 otherwise
%                 opts.globalmin = globalmin (if known)
%                 opts.globalxmin = globalxmin (if known)
%                 opts.dimension = problem dimension
%                 opts.showits   = 1 print iteration status
%                 opts.ep        = global/local weight parameter
%                 opts.tol       = tolerance for termination if
%                                  testflag = 1
%
%       bounds  - (n x 2) matrix of bound constraints LB <= x <= UB
%                 The first column is the LB bounds and the second
%                 column contains the UB bounds
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
% Jones, D.R.: Direct global optimization algorithm. In: Floudas, C.A., 
% Pardalos, P.M. (Eds.) Encyclopedia of Optimization, pp. 431–440. 
% Springer, Boston (2001). https://doi.org/10.1007/0-306-48332-7_93
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[third, VAL, MSS] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS] = Initialization(VAL, OPTI, Problem, MSS);

while VAL.perror > OPTI.TOL                                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    POH = Selection(VAL, MSS, Xmin);                                 %GL

    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL] = Subdivision(VAL, Problem, third, MSS, POH);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS);
end                                                         % End of while

% Return value
minima      = Fmin;
if OPTI.G_nargout == 2
    xatmin    = Xmin;
elseif OPTI.G_nargout == 3
    xatmin    = Xmin;
    history   = VAL.history(1:VAL.itctr, 1:4);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% AUXILIARY FUNCTION BLOCK
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Function  : Selection
% Purpose   : Selection of potential optimal hyper-rectangles
%--------------------------------------------------------------------------
function POH = Selection(VAL, MSS, Xmin)
%--------------------------------------------------------------------------
% Calculate Euclidean Distatnces
Euclid_dist = sum((Xmin(:, 1) - MSS.CC(:, 1:VAL.I)).^2, 1).^0.5;

% Identify potential optimal hyper-rectangles
S = Find_pohas(MSS.FF(1:VAL.I), MSS.DD(1:VAL.I), 1:VAL.I);

D = Find_pohas(Euclid_dist(1:VAL.I), MSS.DD(1:VAL.I), 1:VAL.I);

% Find unique set of potential optimal hyper-rectangles

D(ismember(D, intersect(D, S))) = [];
D = sortrows([D; MSS.DD(D)].', 2).';
S = sortrows([S; MSS.DD(S)].', 2).';
POH = [S(1, :), D(1, :)];
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function:   Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function boxes = Find_pohas(fc, szes, indexs)
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
% Function  : Options
% Purpose   : Get options from inputs
%--------------------------------------------------------------------------
function [OPTI, VAL] = Options(opts, narg, Problem, bounds)
%--------------------------------------------------------------------------
% Determine option values
if nargin < 3 && isempty(opts)
    opts = [];
end
getopts(opts, 'maxits', 1000,'maxevals', 100000, 'maxdeep', 1000,...
    'testflag', 0, 'tol', 0.01, 'showits', 1, 'dimension', 1, 'ep',...
    1e-4, 'globalmin', 0, 'globalxmin', 0);

if isempty(bounds)
% Return the problem information.
    getInfo = feval(Problem.f);
    
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
OPTI.ep        = ep;       % global/local weight parameter
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : BEGIN
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [third, VAL, MSS] = Alocate(OPTI, VAL)
%--------------------------------------------------------------------------
tic                                     % Mesure time
[VAL.time, VAL.Local, VAL.nLs] = deal(0); % initial time
VAL.itctr   = 1;                        % initial iteration
VAL.perror  = 10;                       % initial perror

% alociate MAIN sets
z   = round(OPTI.MAXevals);
MSS = struct('FF', zeros(1, z), 'DD', -ones(1, z),...
    'LL', zeros(VAL.n, z), 'CC', zeros(VAL.n, z));

third       = zeros(1, OPTI.MAXdeep);   % delta values
third(1)    = 1/3;                      % first delta
for i = 2:OPTI.MAXdeep                  % all delta
    third(i)  = (1/3)*third(i - 1);
end
if OPTI.G_nargout == 3
    VAL.history = zeros(OPTI.MAXits, 4); % allocating history
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
VAL.I  = 1;                                        % evaluation counter
MSS.DD(1) = 1;                                     % initial diameter
MSS.CC(:, 1) = ones(VAL.n, 1)/2;                   % initial midpoint
Xmin = abs(VAL.b - VAL.a).*(MSS.CC(:, 1)) + VAL.a; % initial point
[Fmin, MSS.FF(1), VAL.fmin_old, VAL.Flocal] = deal(feval(Problem.f, Xmin));
VAL.minas = zeros(VAL.n, 1);
VAL.x_min_general = Xmin;
Xmin = MSS.CC(:, 1);
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
    VAL.history(VAL.itctr, 1) = 0;
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
[Fmin, ~] =  min(MSS.FF(1:VAL.I));

II         = find(MSS.FF(1:VAL.I) == Fmin);
fminindex  = II(find(MSS.DD(II) == max(MSS.DD(II)), 1, 'last'));
Xmin              = MSS.CC(:, fminindex);


if OPTI.showITS == 1                % Show iteration stats
    VAL.time = toc;
    fprintf(...
    'Iter: %4i  f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, Fmin, VAL.time, VAL.I);
end

if OPTI.TESTflag == 1               % Check for stop condition
    if OPTI.globalMIN ~= 0            % Calculate error if globalmin known
        VAL.perror = 100*(min(VAL.Flocal, Fmin) -...
            OPTI.globalMIN)/abs(OPTI.globalMIN);
    else
        VAL.perror = 100*min(VAL.Flocal, Fmin);
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

if VAL.I + VAL.Local > OPTI.MAXevals % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if max(max(MSS.LL)) + 1 > OPTI.MAXdeep  % Have we exceeded max deep?
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
function [MSS, VAL] = Subdivision(VAL, O, third, MSS, POH)
%--------------------------------------------------------------------------
for i = 1:size(POH, 2)
    % Initial calculations
    ls = find(MSS.LL(:, POH(i)) == min(MSS.LL(:, POH(i))));
    if length(ls) ~= 1
        [~, indexas] = min(VAL.minas(ls));
        ls = ls(indexas);
    end
    VAL.minas(ls) = VAL.minas(ls) + 1;
    delta = third(min(MSS.LL(:, POH(i))) + 1);
    MSS.LL(ls, POH(i)) = MSS.LL(ls, POH(i)) + 1;
    MSS.DD(POH(i)) = 1/2*norm((1/3*(ones(VAL.n, 1))).^(MSS.LL(:, POH(i))));
    
    % Calculate new points and evaluate at the objective function
    [MSS.CC(:, VAL.I + 1), MSS.CC(:, VAL.I + 2)] = deal(MSS.CC(:, POH(i)));
    MSS.CC(ls, VAL.I + 1) = MSS.CC(ls, VAL.I + 1) - delta;
    MSS.CC(ls, VAL.I + 2) = MSS.CC(ls, VAL.I + 2) + delta;

    % Store information
    MSS.FF(VAL.I + 1) = feval(O.f, abs(VAL.b - VAL.a).*MSS.CC(:,...
        VAL.I + 1) + VAL.a);
    MSS.FF(VAL.I + 2) = feval(O.f, abs(VAL.b - VAL.a).*MSS.CC(:,...
        VAL.I + 2) + VAL.a);
    
    [MSS.LL(:, VAL.I + 1), MSS.LL(:, VAL.I + 2)] = deal(MSS.LL(:, POH(i)));
    [MSS.DD(VAL.I + 1), MSS.DD(VAL.I + 2)]       = deal(MSS.DD(POH(i)));
    
    % Update information
    VAL.I = VAL.I + 2;
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