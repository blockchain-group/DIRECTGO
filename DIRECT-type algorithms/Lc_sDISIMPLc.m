function [minima, xatmin, history] = Lc_sDISIMPLc(Problem, bounds, opts)
%--------------------------------------------------------------------------
% Function   : Lc_sDISIMPLc
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 03/17/2020
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = Lc_sDISIMPLc(Problem, bounds, opts)
%       s         - static memory management in data structure
%       DISIMPL   - DISIMPL(dividig simplices) algorithm
%       c         - Eevaluations are performed on the center of the simplices
%       Lc        - Linear constraint handling tehnic
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f       = Objective function handle
%                 Problem.constraint  = constraint functions
%--------------------------------------------------------------------------
%   function [c, ceq] = constraints(x)
%       [c, ceq] = deal([]);
%       c(i)   = g_i(x); inequality constraints i = 1....p
%       ceq(i) = h_i(x); equality constraints i = 1....r   
%   end
%--------------------------------------------------------------------------
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
% Paulavicius, R., and Zilinskas, J. "Simplicial Lipschitz optimization 
% without the Lipschitz constant". 59, 23–40 (2014) Journal of Global 
% Optimization. (2014). DOI 10.1007/s10898-013-0089-3 
%
% Paulavicius, R., Zilinskas, J. Advantages of simplicial partitioning 
% for Lipschitz optimization problems with linear constraints. 
% Optim Lett 10, 237–246 (2016). DOI 10.1007/s11590-014-0772-4
%--------------------------------------------------------------------------

% Get options
OPTI = Options(opts, nargout);

% Alocate sets and create initial variables
[VAL, MSS, V] = Alocate(bounds, OPTI, Problem);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS, Simpl] = Initialization(VAL, OPTI,...
    Problem, MSS, V);

while VAL.perror > OPTI.TOL                                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    POH = Find_po(MSS.FF(1:VAL.m), Fmin, OPTI.ep, MSS.DD(1:VAL.m));

    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL, Simpl] = Subdivision(VAL, Problem, MSS, POH, Simpl);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS, Simpl);
end                                                         % End of while

% Return value
minima      = Fmin;
if OPTI.G_nargout == 2
    xatmin    = (abs(VAL.b - VAL.a)).*Xmin(:, 1) + VAL.a;
elseif OPTI.G_nargout == 3
    xatmin    = (abs(VAL.b - VAL.a)).*Xmin(:, 1) + VAL.a;
    history   = VAL.history(1:(VAL.itctr - 1), 1:4);
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
function OPTI = Options(opts,narg)
%--------------------------------------------------------------------------
% Determine option values
if nargin < 3 && (isempty(opts))
    opts = [];
end
getopts(opts, 'maxits', 1000,'maxevals', 100000, 'maxdeep', 1000,...
    'testflag', 0, 'globalmin', 0, 'tol', 0.01, 'showits', 1,...
    'ep', 1e-4);

OPTI.G_nargout = narg;     % output arguments
OPTI.MAXits    = maxits;   % Fmax of iterations
OPTI.MAXevals  = maxevals; % Fmax # of function evaluations
OPTI.MAXdeep   = maxdeep;  % Fmax number of side divisions
OPTI.TOL       = tol;      % allowable relative error if f_reach is set
OPTI.TESTflag  = testflag; % terminate if global minima is known
OPTI.globalMIN = globalmin;% minimum value of function
OPTI.showITS   = showits;  % print iteration stat
OPTI.ep         = ep;      % global/local weight parameter
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : BEGIN
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [VAL, MSS, V] = Alocate(bounds, OPTI, Problem)
%--------------------------------------------------------------------------
tic                                     % Mesure time
VAL.a       = bounds(:, 1);             % left bound
VAL.b       = bounds(:, 2);             % right bound
VAL.n       = size(bounds, 1);          % dimension
VAL.nn      = size(bounds, 1) + 1;
VAL.alpha   = 1/(VAL.nn);
VAL.time    = 0;                        % initial time
VAL.itctr   = 1;                        % initial iteration
VAL.perror  = 10;                       % initial perror 

V = FindVertices(bounds, Problem);
for i = 1:VAL.n
    VAL.a(i) = min(V(:, i));
    VAL.b(i) = max(V(:, i));
    V(:, i) = (V(:, i)  - VAL.a(i) )/( VAL.b(i) - VAL.a(i));
end

% alociate MAIN sets
z   = round(OPTI.MAXevals);
MSS = struct('FF', zeros(1, z), 'DD', -ones(1, z));

if OPTI.G_nargout == 3
    VAL.history = zeros(OPTI.MAXits, 4); % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [OPTI, VAL, Xmin, Fmin, MSS, Simpl] = Initialization(VAL, OPTI,...
    Problem, MSS, V)
%--------------------------------------------------------------------------
if VAL.n == 2
    dt = delaunayTriangulation(V(:,1),V(:,2)); 
elseif VAL.n == 3
    dt = delaunayTriangulation(V(:,1),V(:,2),V(:,3)); 
else
    dt = delaunayn(V,{'Qt', 'Qbb', 'Qc', 'Qz'});
end

if VAL.n < 4
    VAL.m = size(dt.ConnectivityList, 1);
    for i = 1:VAL.m
        Simpl(i).V = V(dt.ConnectivityList(i, :), :)'; %#ok<AGROW>
    end
else
    VAL.m = size(dt, 1);
    for i = 1:VAL.m
        Simpl(i).V = V(dt(i, :), :)'; %#ok<AGROW>
    end
end

for i = 1:VAL.m
    Simpl(i).SC = sum(VAL.alpha*Simpl(i).V')'; 
    Simpl(i).dSC = 0; 
    for j = 1:VAL.nn
        Simpl(i).dSC = max(Simpl(i).dSC, norm(Simpl(i).SC - Simpl(i).V(:, j)));
    end
    Simpl(i).FV = feval(Problem.f, abs(VAL.b - VAL.a).*Simpl(i).SC + VAL.a);
end  

for i = 1:VAL.m
    % Vector with longest distance from center to vertices
    MSS.DD(i) = Simpl(i).dSC;  
    % Vector with function f(x) values
    MSS.FF(i) = Simpl(i).FV;  
end

% Find vertex with minimal function value
[Fmin, index] = min(MSS.FF(1:VAL.m));  	% initial minima
Xmin = Simpl(index).SC;                 % initial point

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
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
function [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS, Simpl)
%--------------------------------------------------------------------------
[Fmin, fminindex] =  min(MSS.FF(1:VAL.m));
Xmin = Simpl(fminindex).SC;

if OPTI.showITS == 1                % Show iteration stats
    VAL.time = toc;
    fprintf(...
    'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, Fmin, VAL.time, VAL.m);
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

if VAL.m > OPTI.MAXevals       % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr,1) = VAL.itctr;
    VAL.history(VAL.itctr,2) = VAL.m;
    VAL.history(VAL.itctr,3) = Fmin;
    VAL.history(VAL.itctr,4) = VAL.time;
end

VAL.itctr = VAL.itctr + 1 ;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Subdivision
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [MSS, VAL, Simpl] = Subdivision(VAL, O, MSS, POH, Simpl)
%--------------------------------------------------------------------------
for i = 1:size(POH, 2)   
    longDist = 0;
    for i1 = 1:VAL.n
        for j1 = i1 + 1:VAL.nn
            dist = norm(Simpl(POH(i)).V(1:VAL.n, i1) -...
                Simpl(POH(i)).V(1:VAL.n, j1), 2);
            if (dist - longDist > 1e-12)
                longDist = dist;
                kk = i1;
                ll = j1;
            end
        end
    end
    % New vertex point
    midPoint1 = 1/3*Simpl(POH(i)).V(1:VAL.n, kk) +...
        2/3*Simpl(POH(i)).V(1:VAL.n, ll);
    % New vertex point
    midPoint2 = 2/3*Simpl(POH(i)).V(1:VAL.n, kk) +...
        1/3*Simpl(POH(i)).V(1:VAL.n, ll); 
    
    [Simpl(VAL.m + 1).V, Simpl(VAL.m + 2).V] = deal(Simpl(POH(i)).V);
    [Simpl(POH(i)).V(:, ll), Simpl(VAL.m + 2).V(:, kk)] = deal(midPoint1);
    [Simpl(POH(i)).V(:, kk), Simpl(VAL.m + 1).V(:, ll)] = deal(midPoint2); 
    
    Simpl(POH(i)).dSC = roundn(max(sqrt(sum((Simpl(POH(i)).SC -...
        Simpl(POH(i)).V(:, 1:VAL.nn)).^2, 1))), -15);

    % Vector with longest distance from SC to vertices
    MSS.DD(POH(i)) = Simpl(POH(i)).dSC; 
    Simpl(VAL.m + 1).SC = sum(VAL.alpha*Simpl(VAL.m + 1).V')'; 
    
    Simpl(VAL.m + 1).dSC = roundn(max(sqrt(sum((Simpl(VAL.m + 1).SC -...
        Simpl(VAL.m + 1).V(:, 1:VAL.nn)).^2, 1))), -15);

    Simpl(VAL.m + 1).FV = feval(O.f, abs(VAL.b -...
        VAL.a).*Simpl(VAL.m + 1).SC + VAL.a);
    % Vector with longest distance from incenter to vertices
    MSS.DD(VAL.m + 1) = Simpl(VAL.m + 1).dSC;  
    MSS.FF(VAL.m + 1) = Simpl(VAL.m + 1).FV;
    
    Simpl(VAL.m + 2).SC = sum(VAL.alpha*Simpl(VAL.m + 2).V')'; 
    
    Simpl(VAL.m + 2).dSC = roundn(max(sqrt(sum((Simpl(VAL.m + 2).SC -...
        Simpl(VAL.m + 2).V(:, 1:VAL.nn)).^2, 1))), -15);
    
    Simpl(VAL.m + 2).FV = feval(O.f, abs(VAL.b -...
        VAL.a).*Simpl(VAL.m + 2).SC + VAL.a);
           
    MSS.DD(VAL.m + 2) = Simpl(VAL.m + 2).dSC;  
    MSS.FF(VAL.m + 2) = Simpl(VAL.m + 2).FV;
    VAL.m = VAL.m + 2; 
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_po
% Purpose    :  Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function S = Find_po(fc, f_min, epsilon, szes)
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

Su = cell(1, jj);
for i = idx:jj
    idx2 = find((abs(fc - d_min(i)) <= 1E-12) & (szes == d(i)));
    Su{i} = idx2;
end
S_1 = cell2mat(Su);
Su = cell(1, jj);
if length(d)-idx > 1
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
% Function   :  FindVertices
% Purpose    :  Cover hyper-rectangle by simplex
%--------------------------------------------------------------------------
function V = FindVertices(bounds,Problem)
% Find intersecting vertices accorting to linear constrains
warning('off','all');
con_g = feval(Problem.constraint, zeros(1, size(bounds, 1)));
bounds(:, 1) = -bounds(:, 1);
B = [-con_g';reshape(bounds.',1,[])'];
idx = (2*(1:size(bounds, 1))-1)';
Identity = repelem(eye(size(bounds, 1)), 2, 1);

Identity(idx, :) = -Identity(idx, :);
tara = zeros(length(con_g), size(bounds, 1));
%--------------------------------------------------------------------------
for i = 1:size(bounds, 1)
    con = feval(Problem.constraint, Identity(i*2, :));
    tara(:, i) = con - con_g;
end

A = [tara(1:length(con_g), :); Identity];
cn = size(A, 1); % number of constrains
Index = nchoosek(1:cn, size(A, 2));
sz = size(Index, 1);
X = zeros(sz, size(A, 2));
ii = 0;
for i = 1:size(Index, 1)
    XX = roundn(A(Index(i, :), :)\B(Index(i, :)), -12);  
    gg = sum(XX);
    if ~isnan(gg) && gg ~= Inf && gg ~= -Inf
        ii = ii + 1;
        X(ii, :) = XX'; 
    end
end
X = X(1:ii, :);
X = unique(X, 'rows'); % Remove duplicate points

% Verification: Discard solutions outside feasible region
V = []; % Vertex set
for i = 1:size(X, 1)
    corr_cn = true; % Let's say it is inside constrains
    for j = 1:cn 
        if (sum(A(j, :)*X(i, :)') > B(j) + 1e-10)
            corr_cn = false;
            break;
        end
    end
    if corr_cn == true
        V = [V; X(i, :)]; %#ok<AGROW>
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