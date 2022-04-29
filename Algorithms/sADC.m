function [minima, xatmin, history] = sADC(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : sADC
% Written by : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Written by : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 04/29/2020
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = sADC(Problem, opts, bounds)
%       s         - static memory management in data structure
%       ADC       - Adaptive diagonal curves
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
% Sergeyev, Y.D., Kvasov, D.E.: Global search based on efficient diagonal 
% partitions and a set of Lipschitz constants. SIAM J. Optim. 16(3), 
% 910–937 (2006). https://doi.org/10.1137/040621132
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[VAL, MSS, RECT, third] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS, RECT] = Initialization(VAL, OPTI,...
    Problem, MSS, RECT);

while VAL.perror > OPTI.TOL                                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    [POH, VAL] = Find_poh(MSS.FF(1:VAL.I), Fmin, OPTI.ep,...
        MSS.DD(1:VAL.I), VAL, MSS.LL(:, 1:VAL.I));

    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL, RECT, Xmin, Fmin] = Subdivision(VAL, Problem, MSS, POH,...
        RECT, third, Fmin, Xmin);
    
    % Update minima and check stopping conditions
    VAL = Arewedone(OPTI, VAL, MSS, Fmin);
end                                                         % End of while

% Return value
minima      = Fmin;
if OPTI.G_nargout == 2
    xatmin    = (abs(VAL.b - VAL.a)).*Xmin(:, 1) + VAL.a;
elseif OPTI.G_nargout == 3
    xatmin    = (abs(VAL.b - VAL.a)).*Xmin(:, 1) + VAL.a;
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
    1e-5, 'globalmin', 0, 'globalxmin', 0);

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
function [VAL, MSS, RECT, third] = Alocate(OPTI, VAL)
%--------------------------------------------------------------------------
tic                                     % Mesure time
[VAL.time, VAL.global] = deal(0);       % initial time
VAL.itctr   = 1;                        % initial iteration
VAL.perror  = 10;                       % initial perror

% alociate MAIN sets
z   = OPTI.MAXevals;
MSS = struct('FF', nan(1, (z)), 'DD', -ones(1, (z)), 'LL',...
    zeros(VAL.n, (z)), 'CC', zeros(VAL.n, (z)), 'F', zeros(1, (z)));
third       = zeros(1, OPTI.MAXdeep);   % delta values
third(1)    = 1/3;                      % first delta
for i = 2:OPTI.MAXdeep                  % all delta
    third(i)  = (1/3)*third(i - 1);
end

RECT = struct('p1', zeros(VAL.n, 1), 'p2', zeros(VAL.n, 1),...
    'f1', zeros(1, 1), 'f2', zeros(1, 1));

if OPTI.G_nargout == 3
    VAL.history = zeros(OPTI.MAXits, 4); % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [OPTI, VAL, Xmin, Fmin, MSS, RECT] = Initialization(VAL, OPTI,...
    Problem, MSS, RECT)
%--------------------------------------------------------------------------
VAL.I = 1;                                        % evaluation counter
RECT(1).p1 = zeros(VAL.n, 1);                     % initial point p1
RECT(1).p2 = ones(VAL.n, 1);                      % initial point p2 
RECT(1).f1 = feval(Problem.f, VAL.a);
RECT(1).f2 = feval(Problem.f, VAL.b);
MSS.LL(:, 1) = zeros(VAL.n, 1);
MSS(1).CC(:, 1) = zeros(VAL.n, 1);
MSS(1).CC(:, 2) = ones(VAL.n, 1);
VAL.e = 2;
Fmin = min(RECT(1).f1, RECT(1).f2);
    
if RECT(1).f1 < RECT(1).f2
    Xmin = RECT(1).p1;
else
    Xmin = RECT(1).p2;
end
MSS.FF(1) = (RECT(1).f1 + RECT(1).f2)/2;
MSS.DD(1) = 1;
[MSS.F(1), VAL.fMinBeforeImpr] = deal(Fmin);
VAL.fMinNotImpr = 0;

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
    VAL.history(VAL.itctr, 2) = VAL.e;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
function VAL = Arewedone(OPTI, VAL, MSS, Fmin)
%--------------------------------------------------------------------------
if (abs(Fmin - VAL.fMinBeforeImpr) < 0.01*abs(Fmin))
    VAL.fMinNotImpr = VAL.fMinNotImpr + 1;
else
    VAL.fMinNotImpr    = 0;
    VAL.fMinBeforeImpr = Fmin;
end

VAL.time = toc;

if OPTI.showITS == 1                % Show iteration stats
    fprintf(...
    'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, Fmin, VAL.time, VAL.e);
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

if VAL.e > OPTI.MAXevals       % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if max(max(MSS.LL)) > OPTI.MAXdeep  % Have we exceeded max deep?
    disp('Exceeded Max depth. Increse maxdeep');
    VAL.perror = -10;
end

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.e;
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
function [MSS, VAL, R, Xmin, Fmin] = Subdivision(VAL, Problem, MSS,...
    POH, R, third, Fmin, Xmin)
%--------------------------------------------------------------------------
while isempty(POH) == 0
    % Initial calculations
    VAL.I = VAL.I + 2; O = POH(1); A = VAL.I - 1; B = VAL.I;
    
    ls = find(MSS.LL(:, O) == min(MSS.LL(:, O)), 1, 'first');
    delta = third(min(MSS.LL(:, O)) + 1);
    
    MSS.LL(ls, O) = MSS.LL(ls, O) + 1;
    [MSS.LL(:, B), MSS.LL(:, A)] = deal(MSS.LL(:, O));
    
    [MSS.DD(B), MSS.DD(A), MSS.DD(O)] = deal...
        (round(1/2*norm((1/3*(ones(VAL.n, 1))).^(MSS.LL(:, O)), 2), 16));
    
    [R(B), R(A)] = deal(R(O));
    
    % New sampling points
    if R(O).p1(ls) < R(O).p2(ls)
        R(A).p1(ls) = R(A).p1(ls) + 2*delta;  
        R(A).p2(ls) = R(A).p2(ls) - 2*delta;
        R(B).p1 = R(A).p1; R(O).p2 = R(A).p2;
    else
        R(A).p1(ls) = R(A).p1(ls) - 2*delta;  
        R(A).p2(ls) = R(A).p2(ls) + 2*delta; 
        R(B).p1 = R(A).p1; R(O).p2 = R(A).p2; 
    end
    
    idx = find(~any(bsxfun(@minus, MSS.CC(:, 1:VAL.e), R(A).p1)), 1, 'first');
    if isempty(idx)
        VAL.e = VAL.e + 1;
        MSS.CC(:, VAL.e) = R(A).p1;
    end
    idx = find(~any(bsxfun(@minus, MSS.CC(:, 1:VAL.e), R(A).p2)), 1, 'first');
    if isempty(idx)
        VAL.e = VAL.e + 1;
        MSS.CC(:, VAL.e) = R(A).p2;
    end
    % Evaluate function at points
    [R(B).f1, R(A).f1] = deal(feval(Problem.f,...
        (abs(VAL.b - VAL.a).*R(A).p1 + VAL.a)));  
  	[R(A).f2, R(O).f2] = deal(feval(Problem.f,...
        (abs(VAL.b - VAL.a).*R(A).p2 + VAL.a)));  

    % Update min func. value for POH(1) rectangle       
    MSS.FF(O) = (R(O).f1 + R(O).f2)/2; MSS.F(O) = min(R(O).f1, R(O).f2);      
    MSS.FF(A) = (R(A).f1 + R(A).f2)/2; MSS.F(A) = min(R(A).f1, R(A).f2);   
    MSS.FF(B) = (R(B).f1 + R(B).f2)/2; MSS.F(B) = min(R(B).f1, R(B).f2); 
    
    f_calls = [R(O).f1, R(O).f2, R(A).f1, R(A).f2, R(B).f1, R(B).f2];
    pts = [R(O).p1, R(O).p2, R(A).p1, R(A).p2, R(B).p1, R(B).p2];
    [ff, index] = min(f_calls);
    
    if ff < Fmin, Fmin = ff; Xmin = pts(:, index); end
    POH(1) = [];
end
%--------------------------------------------------------------------------
return


%--------------------------------------------------------------------------
% Function   :  Find_po
% Purpose    :  Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function [S, VAL] = Find_poh(fc, f_min, epsilon, szes, VAL, lengths)
% Find all rects on hub
E          = max(epsilon*abs(f_min), 1E-8);
[~, i_min] = min((fc - f_min + E)./szes);
d          = unique(szes);
idx        = find(d == szes(i_min));
jj         = length(d);
d_min      = zeros(1, length(d));

for i = 1:jj
    d_min(i) = min(fc(szes == d(i)));
end

if (VAL.fMinNotImpr > VAL.n + 1)  
    [~, indexas] = min(fc);
    diff_szes = sum(lengths, 1);
    p = diff_szes(indexas);
    q = min(diff_szes); 
    idx = fix((q + p)/2) + 1;
    pp = zeros(VAL.n, jj);
    for i = 1:jj
        pp(:, i) = lengths(:, find(szes == d(i), 1, 'first'));
    end
    idx = find((sum(pp, 1) + 1) == idx);
    VAL.global = VAL.global + 1;
    if VAL.global > 2^(VAL.n + 1)
        [VAL.global, VAL.fMinNotImpr] = deal(0);
    end
end

Su = cell(1, jj);
for i = idx:jj
    idx2 = find((fc - d_min(i) <= 1E-12) & (szes == d(i)));
    Su{i} = idx2;
end
S_1 = cell2mat(Su);
Su = cell(1, jj);
if jj - idx > 1
    a1 = szes(i_min);
    b1 = fc(i_min);
    a2 = d(length(d));
    b2 = d_min(length(d));
    % The line is defined by: y = slope*x + const
    slope = (b2 - b1)/(a2 - a1);
    const = b1 - slope*a1;
    for i = 1 : length(S_1)
        j = S_1(i);
        if fc(j) <= slope*szes(j) + const + 1E-08
            Su{i} = j;
        end
    end
    S_2 = cell2mat(Su);
    if isempty(S_2)
        S_2 = S_1;
    end
    % S_2 now contains all points in S_1 which lies on or below the line
    % Find the points on the convex hull defined by the points in S_2
    xx = szes(S_2);
    yy = fc(S_2);
    h = conhull(xx, yy); % conhull is an internal subfunction
    S_3 = S_2(h);
else
    S_3 = S_1;
end
S = S_3;
return

%--------------------------------------------------------------------------
% Function   :  conhull
% Purpose    :  conhull returns all points on the convex hull.
%--------------------------------------------------------------------------
function h = conhull(x, y)
% conhull returns all points on the convex hull, even redundant ones.
format long;
x = x(:);
y = y(:);
xyAR = roundn([x y], -12);
xyUnique = unique(xyAR, 'rows');
x = xyUnique(:, 1);
y = xyUnique(:, 2);
[w, m, cond] = deal(length(x));
cond = cond - 1;
Z = cell(1, m);
xy = [x, y];
if m == 2
    for i = 1:m
        Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
    end
    h = cell2mat(Z);
    return;
end
if m == 1
    for i = 1:m
        Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
    end
    h = cell2mat(Z);
    return;
end
[v, START, flag] = deal(1);

% Index vector for points in convex hull
h = (1:length(x))';
while (cond ~= START) || (flag == 1)
    if cond == w,  flag = 0; end                                                                    
    a = h(v);                                                                  
    
    if v == m, b = 1; else, b = v + 1; end
    
    b = h(b); 
    if v == m, c = 1; else, c = v + 1; end
    if c == m, c = 1; else, c = c + 1; end
    c = h(c);
    if det([x(a) y(a) 1 ; x(b) y(b) 1 ; x(c) y(c) 1]) >= 0
        leftturn = 1;
    else
        leftturn = 0;
    end
    if leftturn
        if v == m, v = 1; else, v = v + 1; end
        cond = v;
    else
        if v == m, j = 1; else, j = v + 1; end
        h(j:m - 1) = h(j + 1:m);
        m = m - 1;
        w = w - 1;
        if v == 1, v = m; else, v = v - 1; end
        if v == m, cond = 1; else, cond = v + 1; end                        
    end

end
xy = [x(h(1:m)), y(h(1:m))];
hh = size(xy, 1);
Z = cell(1, hh);
for i = 1:hh
    Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
end
h = cell2mat(Z);
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