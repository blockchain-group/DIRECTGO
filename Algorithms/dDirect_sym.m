function [minima, xatmin, history] = dDirect_sym(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : dDirect_sym
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dDirect_sym(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       sym       - DIRECT extension for symmetric functions
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
% Grbic, R., Nyarko, E.K., Scitovski, R. "A modifcation of the DIRECT
% method for Lipschitz global optimization for a symmetric function".
% J. Global Optim. 1–20 (2012). DOI 10.1007/s10898-012-0020-3
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[SS, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[MSS, CE, VAL] = Alocate(SS, VAL);

% Initialization step
[MV, MSS, CE, VAL, minval, xatmin] = Initialization(Problem,...
    MSS, CE, VAL, SS);

while VAL.perror > SS.TOL                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    [POH, VAL] = Find_po(MSS, CE, SS, minval, MV(1, :), VAL);
    
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
    history   = VAL.history(1:VAL.itctr, 1:4);
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
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = Alocate(SS, VAL)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
VAL.time       = 0;                         % initial time
VAL.fcount     = 1;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
CE             = zeros(1, VAL.n*SS.MAXdeep);% collum counter
VAL.m          = 1;
% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
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
    MSS, CE, VAL, SS)
%--------------------------------------------------------------------------
% Create Initial values
MV(1:2,1) 	= 1;                                     % first fake value
MSS(1).L(:, 1) = ones(VAL.n, 1)/2;                   % Lengths
MSS(1).C(:, 1) = ones(VAL.n, 1)/2;                   % Center point
MSS(1).E(1)   = single(1);                           % Index
MSS(1).Diam   = sqrt(sum(MSS(1).L(:, 1).^2));
[MSS(1).F(1), Fmin, MV(1)] = deal(feval(Problem.f, abs(VAL.b -...
    VAL.a).*(MSS(1).C(:, 1)) + VAL.a));
Xmin = MSS(1).C(:, 1);                              % initial point
CE(1) = 1;

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

% Store History
if SS.G_nargout == 3
    VAL.history(VAL.itctr, 1) = 0;
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
                [MSS(i).E, MSS(i).L, MSS(i).F, MSS(i).C] = deal([]);
            end
        else
            C = 1:CE(i); C(POH{i}) = [];
            pp = min(POH{i}):length(C);
            MSS(i).E(pp) = MSS(i).E(C(pp));
            MSS(i).L(:, pp) = MSS(i).L(:, C(pp));
            MSS(i).F(pp) = MSS(i).F(C(pp));
            MSS(i).C(:, pp) = MSS(i).C(:, C(pp));
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
function [MSS, CE, VAL] = STORAGE(Problem, MSS, CE, VAL, PH)
%--------------------------------------------------------------------------
for i = 1:size(MSS, 2)
    if ~isempty(PH{i})
        for j = 1:size(PH{i}, 2)
            [A, TT, VAL, dl] = CALCULS_WORKER(Problem, i, MSS, VAL,...
                PH{i}(j));
            for h = 1:TT
                II              = i + h;               % SET index
                if II > length(CE), CE(II) = 0; end
                MSS(i).L(A.ls(A.order(h)), PH{i}(j)) = dl/2;
                if CE(II) == 0
                    MSS(II).Diam   = sqrt(sum(MSS(i).L(:, PH{i}(j)).^2));
                    [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, 100));
                    [MSS(II).F, MSS(II).E] = deal(zeros(1, 100));
                end
                if CE(II) > size(MSS(II).F, 2)
                    MSS(II).L = [MSS(II).L, zeros(VAL.n, 100)];
                    MSS(II).F = [MSS(II).F, nan(1, 100)];
                    MSS(II).E = [MSS(II).E, zeros(1, 100)];
                    MSS(II).C = [MSS(II).C, zeros(VAL.n, 100)];
                end
                VAL.fcount = VAL.fcount + 2;  
                
                discard = symdirect(VAL, A.c(:, 2*A.order(h) - 1),...
                    MSS(i).L(:, PH{i}(j)));
                [VAL, MSS, CE] = coordi(VAL, discard, MSS, PH{i}(j),...
                    A.c(:, 2*A.order(h) - 1), A.f(2*A.order(h) - 1),...
                    II, CE, i);
                
                discard = symdirect(VAL, A.c(:, 2*A.order(h)),...
                    MSS(i).L(:, PH{i}(j)));
                [VAL, MSS, CE] = coordi(VAL, discard, MSS, PH{i}(j),...
                    A.c(:, 2*A.order(h)), A.f(2*A.order(h)), II, CE, i);
                
                if h == TT
                    discard = symdirect(VAL, MSS(i).C(:, PH{i}(j)),...
                        MSS(i).L(:, PH{i}(j)));
                    if discard == false
                        CE(II)                = CE(II) + 1;
                        MSS(II).F(CE(II))     = A.ff;
                        MSS(II).E(CE(II))     = A.ii;
                        MSS(II).C(:, CE(II))  = A.cc;
                        MSS(II).L(:, CE(II))  = MSS(i).L(:, PH{i}(j));
                    end
                end
            end
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
% Show iteration stats
VAL.time = toc;

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
function [S, VAL] = Find_po(MSS, CE, SS, minval, MV, VAL)
% Find all rects on hub 
  index    = size(MSS, 2);
  d        = [MSS.Diam];
  hulls    = cell(4, index);
  hs       = nan(2, index);
  for i = 1:index
    if CE(i) ~= 0
      ties = find(abs(MSS(i).F(1:CE(i)) -  MV(i)) <= 1E-12);
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

  i_min = find(hs(1, :) == min(hs(1, :))); 
  jj    = find(CE ~= 0, 1, 'first');
  i_m   = i_min - jj;
  jl        = i_m + jj;
  S_1       = cell(1, index);
  S_1(1:jl) = hulls(1, 1:jl); 
  Su        = cell(4, index);
  if i_m > 1        
    a1 = d(i_min);
    b1 = MSS(i_min).F(hs(2, i_min));
    a2 = d(jj);
    b2 = MV(jj);
% The line is defined by: y = slope*x + const
    slope = (b2 - b1)/(a2 - a1);
    const = b1 - slope*a1;
    for i = 1:jl
      if CE(i) ~= 0
        TT = find(MSS(i).F(hulls{1, i}) <= slope*d(i) + const + 1e-12);
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
      if CE(i) ~= 0 
        S{i} = invers(1, invers(4, :) == i);
      end
    end
  else
    S = S_1;
  end
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
    if cond == w,  flag = 0; end                                            % stoping condition                            
    a = h(v);                                                                  % nuo galo einantis elementas
    
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
        if v == m, cond = 1; else, cond = v + 1; end                        % stoping condition
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
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, LL, VAL, dl] = CALCULS_WORKER(Const, i, MSS, VAL, PH)
%--------------------------------------------------------------------------
% Create and allocating empty sets
A.cc      = MSS(i).C(:, PH);
A.ii      = MSS(i).E(PH);
A.ff      = MSS(i).F(PH);

max_L  = max(MSS(i).L(:, PH));
A.ls   = find(MSS(i).L(:, PH) == max_L);
LL     = length(A.ls);
dl     = 2*max_L/3;
A.c    = MSS(i).C(:, PH)*ones(1, 2*LL);

A.c(A.ls, 1:2:end) = A.c(A.ls, 1:2:end) + diag(repelem(dl, LL));
A.c(A.ls, 2:2:end) = A.c(A.ls, 2:2:end) - diag(repelem(dl, LL));
point = abs(VAL.b - VAL.a).*A.c + VAL.a;
A.f = arrayfun(@(x) feval(Const.f, point(:, x)), (1:2*LL));
[~, A.order] = sort([min(A.f(1:2:end), A.f(2:2:end))' A.ls], 1);
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
% Function  : symdirect
% Purpose   : symDIRECT strategy
%--------------------------------------------------------------------------
function discard = symdirect(VAL, c_tp, l_tp)
%--------------------------------------------------------------------------
% Generate all rectangle vertices in terms of O and 1
ok = 1; % true
b = zeros(VAL.n, 1);
k = 1;
V(:, k) = b;
k = k + 1;
while (ok == 1)
    j = VAL.n;
    while ((j > 0) && (b(j) == 1))
        b(j) = 0;
        j = j - 1;
    end
    if (j > 0)
        b(j) = 1;
        V(:, k) = b';
        k = k + 1;
    else
        ok = 0;
    end
end

% Transform to real coordinates
for i = 1:2^VAL.n
    for j = 1:VAL.n
        if V(j, i) == 0
            V(j, i) = c_tp(j, 1) - l_tp(j, 1);
        else
            V(j, i) = c_tp(j, 1) + l_tp(j, 1);
        end
    end
end

% Checking the symmetry for first rectangle
% Avoid distances between two cluster centers smaller than delta
discard = true;
for i = 1:2^VAL.n
    check = 0;
    for j = 1:VAL.n - 1
        if V(j ,i) >= V(j + 1, i); check = check + 1; end
    end
    if check == VAL.n - 1; discard = false; break; end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : coordi
% Purpose   : symDIRECT strategy
%--------------------------------------------------------------------------
function [VAL, MSS, CE] = coordi(VAL, discard, MSS, index, cc,...
    ff, col, CE, i)
%--------------------------------------------------------------------------
if discard == false
    VAL.m                = VAL.m + 1;
    CE(col)              = CE(col) + 1;
    IS                   = CE(col);            % Left index
    MSS(col).L(:, IS)    = MSS(i).L(:, index);
    MSS(col).E(IS)       = VAL.m;
    MSS(col).C(:, IS)    = cc;
    MSS(col).F(IS)       = ff;
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