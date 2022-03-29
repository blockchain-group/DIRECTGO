function [minima, xatmin, history] = Lc_dDISIMPLv(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : Lc_dDISIMPLv
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = Lc_dDISIMPLv(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
%       DISIMPL   - DISIMPL(dividig simplices) algorithm
%       v         - Eevaluations are performed on the vertices of the simplices
%       Lc        - Linear constraint handling tehnic
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f       = Objective function handle
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
%                 opts.ep        = global/local weight parameter
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
% Paulavicius, R., and Zilinskas, J. "Simplicial Lipschitz optimization 
% without the Lipschitz constant". 59, 23–40 (2014) Journal of Global 
% Optimization. (2014). DOI 10.1007/s10898-013-0089-3 
%
% Paulavicius, R., Zilinskas, J. Advantages of simplicial partitioning 
% for Lipschitz optimization problems with linear constraints. 
% Optim Lett 10, 237–246 (2016). DOI 10.1007/s11590-014-0772-4
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[SS, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[MSS, CE, VAL, V] = Alocate(SS, VAL, Problem);

% Initialization step
[MV, MSS, CE, VAL, minval, xatmin] = Initialization(Problem,...
    MSS, CE, VAL, SS, V);

while VAL.perror > SS.TOL                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    POH = Find_po(MSS, CE, SS, minval, MV(1, :));

    % Subdivide potential optimalhyper-rectangles
    [VAL, MSS, CE, MV, minval, xatmin] = Calulcs(VAL, Problem,...
        MSS, CE, POH);

    % Update minima and check stopping conditions
    [VAL, SS] = Arewedone(minval, VAL, SS);
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
    1e-8, 'globalmin', 0, 'globalxmin', 0, 'ep', 1e-4);

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
OPTI.ep        = ep;       % global/local weight parameter
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL, V] = Alocate(SS, VAL, Problem)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
VAL.nn         = VAL.n + 1;
VAL.time       = 0;                         % initial time
VAL.fcount     = 0;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
CE             = zeros(1, VAL.n*SS.MAXdeep);% collum counter

V = FindVertices([VAL.a, VAL.b], Problem);
for i = 1:VAL.n
    VAL.a(i) = min(V(:, i));
    VAL.b(i) = max(V(:, i));
    V(:, i) = (V(:, i)  - VAL.a(i) )/( VAL.b(i) - VAL.a(i));
end

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.nn, 1),...
    'L', zeros(VAL.n, 1));

if SS.G_nargout == 3
    VAL.history = zeros(SS.MAXits, 4);       % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [MV, MSS, CE, VAL, Fmin, Xmin] = Initialization(Problem,...
    MSS, CE, VAL, SS, V)
%--------------------------------------------------------------------------
% Create Initial values
MV(1:2, 1) = 1;                                     % first fake value
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
        MSS(1).Simpl(i).V = V(dt.ConnectivityList(i, :), :)'; 
    end
else
    VAL.m = size(dt, 1);
    for i = 1:VAL.m
        MSS(1).Simpl(i).V = V(dt(i, :), :)'; 
    end
end
V = V';
CE(1) = VAL.m;

% Calculate function values on all vertices V
for i = 1:length(V)     
    V(VAL.nn, i) = feval(Problem.f, abs(VAL.b - VAL.a).*(V(1:VAL.n, i)) + VAL.a);
    VAL.fcount  = VAL.fcount + 1;
end

if VAL.m == 1 
    VAL.fcount = VAL.nn; 
end

% assign function values from the vertex set V
for j = 1:VAL.m 
    for k = 1:VAL.nn 
        for i = 1:length(V)
            if (norm(MSS(1).Simpl(j).V(1:VAL.n, k) - V(1:VAL.n, i),2) == 0)
                MSS(1).Simpl(j).V(VAL.nn, k) = V(VAL.nn, i);
                break;
            end
        end
    end
end
% Find vertex with minimal function value
for i = 1:VAL.m
    fmin = min(MSS(1).Simpl(i).V(VAL.nn, :));  
    MSS(1).F(i) = fmin;     
end

MSS(1).Diam = roundn(sqrt(sum(ones(VAL.n, 1).^2)), -12);
MSS(1).C(:, 1:length(V)) = V;
[Fmin, fminindex] = min(MSS.C(VAL.nn, 1:length(V)));   % initial minima
Xmin = MSS.C(1:VAL.n , fminindex);                            % initial point  

% Check stop condition if global minima is known
if SS.TESTflag  == 1
    if SS.globalMIN ~= 0
        VAL.perror = 100*(Fmin - SS.globalMIN)/abs(SS.globalMIN);
    else
        VAL.perror = 100*Fmin;
    end
else
    VAL.perror   = 2;
end

if SS.G_nargout == 3
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.fcount;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [VAL, MSS, CE, MV, minval, xatmin] = Calulcs(VAL,...
    Problem, MSS, CE, POH)
%--------------------------------------------------------------------------
[MSS, CE, VAL] = STORAGE(Problem, MSS, CE, VAL, POH);

for i = size(POH, 2):-1:1                
    if ~isempty(POH{i})         
        if (CE(i) - size(POH{i}, 2)) == 0
            if find(CE ~= 0, 1, 'first') == i
                [MSS(i).Simpl, MSS(i).F] = deal([]);
            end
        else
            C = 1:CE(i); C(POH{i}) = [];
            pp = min(POH{i}):length(C);
            MSS(i).Simpl(pp) = MSS(i).Simpl(C(pp));
            MSS(i).F(pp) = MSS(i).F(C(pp));
        end
        CE(i) = CE(i) - size(POH{i}, 2);
    end
end

I = 1:size(MSS, 2);
[~, II] = sort([MSS.Diam], 'descend');
MSS(I) = MSS(II);
CE(I) = CE(II);

% Find minima values
[MV, minval, xatmin] = Find_min(CE, MSS, VAL);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = STORAGE(Problem, MSS, CE, VAL, PH)
%--------------------------------------------------------------------------
diam = [MSS.Diam];
for i = 1:size(MSS, 2)
    if ~isempty(PH{i})
        for ii = 1:size(PH{i}, 2)
            % Increase number of simplices
            [Simpl(1), Simpl(2)] = deal(MSS(i).Simpl(PH{i}(ii)));            
            VAL.m = VAL.m + 1;
            
            [longDist, kk, ll, longDist1, longDist2] = deal(0);
            for j = 1:VAL.n
                [dist, jj] = max(sqrt(sum((MSS(i).Simpl(PH{i}(ii)).V(1:VAL.n, j)...
                    - MSS(i).Simpl(PH{i}(ii)).V(1:VAL.n, (j + 1):VAL.nn)).^2, 1)));
                 if dist > longDist
                     longDist = dist;
                     kk = j;
                     ll = jj + j;
                 end
            end
                       
            % Midpoint
            midPoint = (MSS(i).Simpl(PH{i}(ii)).V(1:VAL.n, kk) +...
                MSS(i).Simpl(PH{i}(ii)).V(1:VAL.n, ll))/2;

            idx = find(~any(bsxfun(@minus, MSS(1).C(1:VAL.n,...
                1:VAL.fcount), midPoint)), 1, 'first');
            
            % Calculate function value in new point and add this point to vertex set
            if isempty(idx)
                midPoint(VAL.nn, 1) = feval(Problem.f, abs(VAL.b -...
                    VAL.a).*midPoint + VAL.a);
                VAL.fcount = VAL.fcount + 1;
                MSS(1).C(:, VAL.fcount) = midPoint;
            else
                midPoint(VAL.nn, 1) = MSS(1).C(VAL.nn, idx);
            end
            Simpl(1).V(:, kk) = midPoint;
            Simpl(2).V(:, ll) = midPoint;          
            for j = 1:VAL.n
                dist1 = roundn(max(sqrt(sum((Simpl(1).V(1:VAL.n, j)...
                    - Simpl(1).V(1:VAL.n, (j + 1):VAL.nn)).^2, 1))), -12);
                if dist1 > longDist1
                    longDist1 = dist1;
                end
                dist2 = roundn(max(sqrt(sum((Simpl(2).V(1:VAL.n, j)...
                    - Simpl(2).V(1:VAL.n, (j + 1):VAL.nn)).^2, 1))), -12);
                if dist2 > longDist2
                    longDist2 = dist2;
                end
            end
            
            I = find(diam(1:end) == longDist1, 1, 'first');
            if I > length(CE), CE(I) = 0; end
            if ~isempty(I)
                IC = CE(I) + 1;
                CE(I) = IC;
            else
                diam = [diam, longDist1]; %#ok<AGROW>
                I = length(diam);
                if I > length(CE), CE(I) = 0; end
                IC = CE(I) + 1;
                CE(I) = IC;
                MSS(I).Diam = longDist1;
            end
            MSS(I).Simpl(IC) = Simpl(1);
            MSS(I).F(IC) = min(MSS(I).Simpl(IC).V(VAL.nn, :)); 
            
            I = find(diam(1:end) == longDist2, 1, 'first');
            if I > length(CE), CE(I) = 0; end
            if ~isempty(I)
                IC = CE(I) + 1;
                CE(I) = IC;
            else
                diam = [diam, longDist2]; %#ok<AGROW>
                I = length(diam);
                if I > length(CE), CE(I) = 0; end
                IC = CE(I) + 1;
                CE(I) = IC;
                MSS(I).Diam = longDist2;
            end
            MSS(I).Simpl(IC) = Simpl(2);
            MSS(I).F(IC) = min(MSS(I).Simpl(IC).V(VAL.nn, :));
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
%--------------------------------------------------------------------------
function [VAL, SS] = Arewedone(minval, VAL, SS)
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
% Function  : Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function S = Find_po(MSS, CE, SS, minval, MV)
% Find all rects on hub
index    = size(MSS, 2);
d        = [MSS.Diam];
hulls    = cell(3, index);
hs       = nan(2, index);

for i = 1:index
    if CE(i) ~= 0
        ties     = find((MSS(i).F(1:CE(i)) - MV(i)) <= 1E-12);
        if ~isempty(ties)
            hulls{1, i} = ties;
            hulls{2, i} = d(i)*ones(1, size(ties, 2));
            hulls{3, i} = MSS(i).F(ties);
            hulls{4, i} = i*ones(1, size(ties, 2));
            [hs(1, i), hs(2, i)] = min((MSS(i).F(1:CE(i)) - minval +...
                max(SS.ep*abs(minval), 1E-8))/d(i));
        end
    end
end

jj    = find(~cellfun(@isempty, hulls(1, :)), 1, 'first');
i_min = find(hs(1, :) == min(hs(1, :)));
d_min = MV(1, :);
S_1   = cell(1, index);
S_1(1:i_min)   = hulls(1, 1:i_min);
Su    = cell(4, index);

if i_min - jj > 1
    a1 = d(i_min);
    b1 = MSS(i_min).F(hs(2, i_min));
    a2 = d(jj);
    b2 = d_min(jj);
    % The line is defined by: y = slope*x + const
    slope = (b2 - b1)/(a2 - a1);
    const = b1 - slope*a1;
    for i = 1:i_min
        if ~isempty(hulls{3, i})
            TT = find(hulls{3, i} <= slope*d(i) + const + 1E-08);
            if ~isempty(TT)
                Su{1, i} = hulls{1, i}(TT);
                Su{2, i} = hulls{2, i}(TT);
                Su{3, i} = hulls{3, i}(TT);
                Su{4, i} = i*ones(1, size(TT, 2));
            end
        end
    end
    
    S_2 = cell2mat(Su(1:4, :));
    if isempty(S_2)
        S_2 = cell2mat(hulls(1:4, :));
        S_2 = sortrows(S_2.',2).';
    else
        S_2 = sortrows(S_2.',2).';
    end
    
    % S_2 now contains all points in S_1 which lies on or below the line
    % Find the points on the convex hull defined by the points in S_2
    xx = S_2(2, :);
    yy = S_2(3, :);
    h = conhull(xx, yy); % conhull is an internal subfunction
    invers = S_2(:, h);
    S = cell(1, index);
    for i = 1:index
        if ~isempty(hulls{3, i})
            S{i} = invers(1, invers(4, :) == i);
        end
    end
else
    S = S_1;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_min
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function [MV, minval, xatmin] = Find_min(CE, MSS, VAL)
%--------------------------------------------------------------------------
MV = nan(2, size(MSS, 2));

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        [MV(1, i), MV(2, i)] = min(MSS(i).F(1:CE(i)));
    end
end

[minval, Least] = min(MV(1, :));
xatmin = MSS(Least).Simpl(MV(2, Least)).V(1:VAL.n, 1); 
%--------------------------------------------------------------------------
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