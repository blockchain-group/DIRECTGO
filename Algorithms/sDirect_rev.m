function [minima, xatmin, history] = sDirect_rev(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : sDirect_rev
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 03/17/2020
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = sDirect_rev(Problem, opts, bounds)
%       s         - static memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       rev       - enhancements to reduce the global drag
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
    POH = Find_po(MSS.FF(1:VAL.I), MSS.LL(:,1:VAL.I), Fmin, OPTI.ep,...
        MSS.DD(1:VAL.I));
    
    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL] = Subdivision(VAL, Problem, third, MSS, POH);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin] = Arewedone(Problem, OPTI, VAL, MSS);
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
function [VAL, Fmin, Xmin] = Arewedone(Problem, OPTI, VAL, MSS)
%--------------------------------------------------------------------------
[Fmin, fminindex] =  min(MSS.FF(1:VAL.I));
Xmin              = abs(VAL.b - VAL.a).*MSS.CC(:, fminindex) + VAL.a;

if (abs(Fmin - VAL.fmin_old) > 0.0001) && Fmin < VAL.Flocal     
    options = optimoptions('fmincon', 'StepTolerance', 10^(-10),...
        'Display', 'none');
    [Xmin_local, Fmin_local, ~, output] = fmincon(Problem.f, Xmin, [], [],...
        [], [], VAL.a, VAL.b, [], options);
    VAL.nLs = VAL.nLs + 1;
    VAL.Local = VAL.Local + output.funcCount;
    if (Fmin_local < Fmin)
        [Fmin, VAL.Flocal] = deal(Fmin_local);
        VAL.x_min_general = Xmin_local;
        Xmin = Xmin_local;
    end
end
[Fmin, indexs] = min([Fmin, VAL.Flocal]);
XX           = [Xmin, VAL.x_min_general];
Xmin         = XX(:, indexs);

VAL.fmin_old = Fmin;

if OPTI.showITS == 1                % Show iteration stats
    VAL.time = toc;
    fprintf(...
    'Iter: %4i  f_min: %15.10f    time(s): %10.05f  local searches: %4i    fn evals: %8i\n',...
        VAL.itctr, min(VAL.Flocal, Fmin), VAL.time, VAL.nLs, VAL.I + VAL.Local);
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
    VAL.history(VAL.itctr, 2) = VAL.I + VAL.Local;
    VAL.history(VAL.itctr, 3) = min(VAL.Flocal, Fmin);
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
% Function   :  Find_po                                                                              
% Purpose    :  Return list of PO hyperrectangles                 
%--------------------------------------------------------------------------
function final_pos = Find_po(fc, lengths, minval, ep, szes)
% Find all rects on hub
diff_szes = sum(lengths, 1);
tmp_max = max(diff_szes);
j = 0;
hull = zeros(1, size(fc, 2));
for i = 1:tmp_max + 1
    tmp_idx = find(diff_szes == i - 1);
    iii = min(fc(tmp_idx));
    hullidx = find(fc(tmp_idx) == iii, 1, 'last');
    if ~isempty(hullidx)
        j = j + 1;
        hull(j) = tmp_idx(hullidx);
    end
end
hull = hull(1:j);

% Compute lb and ub for rects on hub
lbound = calc_lbound(lengths, fc, hull, szes);
ubound = calc_ubound(lengths, fc, hull, szes);

% Find indeces of hull who satisfy
maybe_po = find(lbound - ubound <= 0);

if minval ~= 0
    po = find((minval-fc(hull(maybe_po)))./abs(minval) +...
        szes(hull(maybe_po)).*ubound(maybe_po)./abs(minval) >= ep);
else
    po = find(fc(hull(maybe_po)) -...
        szes(hull(maybe_po)).*ubound(maybe_po) <= 0);
end
if isempty(po)
    [~, po] = max(szes(hull(maybe_po)));
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