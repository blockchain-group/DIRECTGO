function [minima, xatmin, history] = sGb_BIRECT(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : sGb_BIRECT
% Written by : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Written by : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 04/29/2020
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = sGb_BIRECT(Problem, opts, bounds)
%       s         - static memory management in data structure
%       BIRECT    - BIRECT(BIsecting RECTangles) algorithm
%       Gb        - globally-biased
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
%                 opts.tighterLB = Use more expensive tighter LB (true/false)
%                 opts.tol       = tolerance for termination if
%                                  testflag = 1
%                 opts.balancedSampling = Perform balanced sampling (true/false)
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
% Paulavicius, R., Chiter, L., and Zilinskas, J. "Global optimization 
% based on bisection of rectangles, function values at diagonals, and a
% set of Lipschitz constants" Journal of Global Optimization. (2018). 
% DOI 10.1007/s10898-016-0485-6 
%
% Paulavicius, R., Sergeyev, Y. D., Kvasov, D. E., & Zilinskas, J. 
% "Globally-biased BIRECT algorithm with local accelerators for expensive 
% global optimization" Expert Systems with Applications. (2020).
% DOI 10.1016/j.eswa.2019.113052
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[VAL, MSS, RECT] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS, RECT] = Initialization(VAL, OPTI,...
    Problem, MSS, RECT);

while VAL.perror > OPTI.TOL                                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    POH = Find_po(MSS.FF(1:VAL.I), Fmin, OPTI.ep, MSS.DD(1:VAL.I), VAL);
    
    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL, RECT] = Subdivision(VAL, Problem, MSS, POH, RECT, OPTI);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS, RECT);
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
    'testflag', 0, 'tol', 0.01, 'showits', 1, 'dimension', 1, 'ep', 1e-4,...
    'balancedSampling', false, 'tighterLB', false, 'globalmin', 0,...
    'globalxmin', 0);

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
OPTI.balancedSampling = balancedSampling; % Perform balanced sampling
OPTI.tighterLB = tighterLB;% Use more expensive tighter LB
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : BEGIN
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [VAL, MSS, RECT] = Alocate(OPTI, VAL)
%--------------------------------------------------------------------------
tic                                     % Mesure time
VAL.time    = 0;                        % initial time
VAL.itctr   = 1;                        % initial iteration
VAL.perror  = 10;                       % initial perror

% alociate MAIN sets
z   = round(OPTI.MAXevals);
MSS = struct('FF', zeros(1, (z/2)), 'DD', -ones(1, (z/2)), 'LL',...
    zeros(VAL.n, (z/2)));

RECT = struct('p1', zeros(VAL.n, 1), 'p2', zeros(VAL.n, 1),...
    'c', zeros(VAL.n, 1), 'f1', zeros(1, 1), 'f2',...
    zeros(1, 1), 'fmin', zeros(1, 1), 'pmin', zeros(VAL.n, 1));

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
RECT(1).p1 = ones(VAL.n, 1)./3;                   % initial point p1
RECT(1).p2 = ones(VAL.n, 1).*(2/3);               % initial point p2
RECT(1).c = (RECT(1).p1 + RECT(1).p2)/2;         
RECT(1).f1 = feval(Problem.f, abs(VAL.b - VAL.a).*RECT(1).p1 + VAL.a);
RECT(1).f2 = feval(Problem.f, abs(VAL.b - VAL.a).*RECT(1).p2 + VAL.a);

if RECT(1).f1 < RECT(1).f2
    RECT(1).fmin = RECT(1).f1;
    RECT(1).pmin = RECT(1).p1;
else
    RECT(1).fmin = RECT(1).f2;
    RECT(1).pmin = RECT(1).p2;
end
Xmin = RECT(1).pmin;

[MSS.FF(1), Fmin, VAL.fMinBeforeImpr] = deal(RECT(1).fmin);
MSS.LL(:, 1) = ones(VAL.n, 1);

VAL.fMinNotImpr = 0;

if OPTI.tighterLB == true 
    % Calculate more expensive but tighter diameter
    RECT(1).v = VertexGeneration(VAL.n); 
    MSS.DD(1) = Diameter(RECT(1).v, RECT(1).p1, RECT(1).p2);
else
    MSS.DD(1) = roundn((2/3)*sqrt(sum(MSS.LL(:, 1).^2)), -12);
end

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
    VAL.history(VAL.itctr, 2) = VAL.I*2;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
function [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS, RECT)
%--------------------------------------------------------------------------
[Fmin, fminindex] =  min(MSS.FF(1:VAL.I));
Xmin              = RECT(fminindex).pmin;

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
        VAL.itctr, Fmin, VAL.time, VAL.I*2);
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

if VAL.I*2 > OPTI.MAXevals       % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if max(max(MSS.LL)) > OPTI.MAXdeep  % Have we exceeded max deep?
    disp('Exceeded Max depth. Increse maxdeep');
    VAL.perror = -10;
end

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.I*2;
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
function [MSS, VAL, RECT] = Subdivision(VAL, O, MSS, POH, RECT, OPTI)
%--------------------------------------------------------------------------
while isempty(POH) == 0
    % Initial calculations
    
    VAL.I = VAL.I + 1;
    max_L = max(MSS.LL(:, POH(1)));
    lsa = find(MSS.LL(:, POH(1)) == max_L);
    ls = lsa(1);
    dl = max_L/2;
    MSS.LL(:, VAL.I) = MSS.LL(:, POH(1));
    [MSS.LL(ls, VAL.I), MSS.LL(ls, POH(1))] = deal(dl);
    
    RECT(VAL.I) = RECT(POH(1));
    % New sampling points
    if RECT(POH(1)).p1(ls) < RECT(POH(1)).p2(ls)
        RECT(VAL.I).p1(ls) = RECT(POH(1)).p1(ls) + dl;  
        RECT(POH(1)).p2(ls) = RECT(VAL.I).p2(ls) - dl; 
    else
        RECT(VAL.I).p1(ls) = RECT(POH(1)).p1(ls) - dl;  
        RECT(POH(1)).p2(ls) = RECT(VAL.I).p2(ls) + dl;  
    end
    
    % Update center points for both rectangles
    RECT(POH(1)).c = (RECT(POH(1)).p1 + RECT(POH(1)).p2)/2;
    RECT(VAL.I).c = (RECT(VAL.I).p1 + RECT(VAL.I).p2)/2;
    
    % Evaluate function at points
    RECT(VAL.I).f1 = feval(O.f, abs(VAL.b - VAL.a).*RECT(VAL.I).p1 + VAL.a);  
  	RECT(POH(1)).f2 = feval(O.f, abs(VAL.b - VAL.a).*RECT(POH(1)).p2 + VAL.a);  

    % Find/update the minimal obj. value and the minimum point
    if RECT(POH(1)).f1 < RECT(POH(1)).f2
        RECT(POH(1)).fmin = RECT(POH(1)).f1;
        RECT(POH(1)).pmin = RECT(POH(1)).p1;
    else
        RECT(POH(1)).fmin = RECT(POH(1)).f2;
        RECT(POH(1)).pmin = RECT(POH(1)).p2;
    end
    
    if RECT(VAL.I).f1 < RECT(VAL.I).f2
        RECT(VAL.I).fmin = RECT(VAL.I).f1;
        RECT(VAL.I).pmin = RECT(VAL.I).p1;
    else
        RECT(VAL.I).fmin = RECT(VAL.I).f2;
        RECT(VAL.I).pmin = RECT(VAL.I).p2;
    end
    
    % Update min func. value for POH(1) rectangle       
    MSS.FF(POH(1)) = RECT(POH(1)).fmin;     
    % Add new min value of the new rectangle
    MSS.FF(VAL.I) = RECT(VAL.I).fmin;  
    
    if OPTI.tighterLB == true
        if (RECT(POH(1)).c(ls) < RECT(VAL.I).c(ls))
            [RECT(POH(1)).v, RECT(VAL.I).v] = UpdateVertices...
                (RECT(POH(1)).v, ls);
        else
            [RECT(VAL.I).v, RECT(POH(1)).v] = UpdateVertices...
                (RECT(POH(1)).v, ls);
        end
        % Calculate diameter for the new rectangles
        [MSS.DD(VAL.I), MSS.DD(POH(1))] = deal(roundn(Diameter...
            (RECT(POH(1)).v, RECT(POH(1)).p1, RECT(POH(1)).p2), -12));
    else
    [MSS.DD(VAL.I), MSS.DD(POH(1))] = deal(roundn((2/3)*sqrt...
        (sum(MSS.LL(:, POH(1)).^2)), -12));
    end
    if (OPTI.balancedSampling == true) && (length(lsa) > 1)
        POH = union(VAL.I, POH);
    else
        POH(1) = [];
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_po
% Purpose    :  Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function S = Find_po(fc, f_min, epsilon, szes, VAL)
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

if (VAL.fMinNotImpr > 2*VAL.n) && mod(VAL.fMinNotImpr, 10) ~= 0
    idx = idx + floor(0.75*(length(d) - idx));            
end

Su = cell(1, jj);
for i = idx:jj
    idx2 = find((fc - d_min(i) <= 1E-12) & (szes == d(i)));
    Su{i} = idx2;
end
S_1 = cell2mat(Su);
Su = cell(1, jj);
if length(d) - idx > 1
    a1 = szes(i_min);
    b1 = fc(i_min);
    a2 = d(length(d));
    b2 = d_min(length(d));
    % The line is defined by: y = slope*x + const
    slope = (b2-b1)/(a2-a1);
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
% Function   :  VertexGeneration
% Purpose    :  Generate all vertices of the unit hyper-rectangle (cube)
%--------------------------------------------------------------------------
function V = VertexGeneration(n)
ok = 1;  % true
b = zeros(1, n);
k = 1;
V(:, k) = b';
k = k + 1;
while (ok == 1)
    j = n;
    while ((j > 0) && (b(j) == 1))
        b(j) = 0;
        j = j - 1;
    end
    if (j > 0)
        b(j) = 1;
        V(:, k) = b';
        k = k + 1;
    else
        ok = 0;  % false
    end
end
return

%--------------------------------------------------------------------------
% Function   :  UpdateVertices
% Purpose    :  Return updated (after branching) coordinates for the 
%               vertices. Branching is done on the idx variable
%--------------------------------------------------------------------------
function [V1, V2] = UpdateVertices(V, idx)
row_max = max(V(idx, :));
row_min = min(V(idx, :));

% Finds maximum/minimum indices of coresponding vertices
I1 = V(idx, :) == row_max;
I2 = V(idx, :) == row_min;

V1 = V; 
V2 = V1;
V1(idx, I1) = (row_min + row_max)/2;  % Shrink from right to left
V2(idx, I2) = (row_min + row_max)/2;  % Shrink from left to right
return

%--------------------------------------------------------------------------
% Function   :  Diameter
% Purpose    :  Calculate diameter for a certain hyper-rectangle
%--------------------------------------------------------------------------
function D = Diameter(V, p1, p2)
% Iterate through all the rectangle vertices
D = 0;  % Set initial diameter's value
for i = 1:size(V, 2)
    D = max(D, min(norm(V(:, i) - p1,2), norm(V(:, i) - p2, 2)));
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