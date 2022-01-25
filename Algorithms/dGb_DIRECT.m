function [minima, xatmin, history] = dGb_DIRECT(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : dGb_DIRECT
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dGb_DIRECT(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
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
% Paulavicius, R., Sergeyev, Y.D., Kvasov, D.E., Zilinskas, J.
% "Globally-biased DISIMPL algorithm for expensive global optimization".
% Journal of Global Optimization 59, 545–567. (2014)
% DOI 10.1007/s10898-014-0180-4.
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
VAL.fMinNotImpr = 0;
VAL.Yaro_epsilon = 0.01;
CE             = zeros(1, VAL.n*SS.MAXdeep);      % collum counter

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
MV(1:2, 1) = 1;                                     % first fake value
MSS(1).L(:, 1) = ones(VAL.n, 1)/2;                  % Lengths
MSS(1).Diam   = roundn(sqrt(sum(MSS.L(:, 1).^2)), -15);
MSS(1).C(:, 1) = ones(VAL.n, 1)/2;               	% Center point
MSS(1).E(1) = single(1);                            % Index
[MSS(1).F(1), Fmin, MV(1), VAL.fMinBeforeImpr] = deal(feval(Problem.f,...
    abs(VAL.b - VAL.a).*(MSS(1).C(:, 1)) + VAL.a));
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
            [DD, TT, VAL] = CALCULS_WORKER(Problem, i, MSS, VAL, PH{i}(j));
            for h = 1:TT
                II = i + h;
                if II > length(CE), CE(II) = 0; end
                if CE(II) == 0
                    MSS(II).Diam = roundn(sqrt(sum(DD.lL(:, h).^2)), -15);
                    [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, 100));
                    [MSS(II).F, MSS(II).E] = deal(zeros(1, 100));
                end
                IL = CE(II) + 1;
                IR = CE(II) + 2;
                CE(II) = IR;
                if CE(II) > size(MSS(II).F, 2)
                    MSS(II).L = [MSS(II).L, zeros(VAL.n, 100)];
                    MSS(II).F = [MSS(II).F, nan(1, 100)];
                    MSS(II).E = [MSS(II).E, zeros(1, 100)];
                    MSS(II).C = [MSS(II).C, zeros(VAL.n, 100)];
                end
                MSS(II).F(IL) = DD.L(DD.lsx(h, 1));
                MSS(II).F(IR) = DD.R(DD.lsx(h, 1));
                MSS(II).E(IL) = DD.eL(h);
                MSS(II).E(IR) = DD.eR(h);
                MSS(II).L(:, IL) = DD.lL(:, h);
                MSS(II).L(:, IR) = DD.lL(:, h);
                MSS(II).C(:, IL) = DD.cL(:, DD.lsx(h, 1));
                MSS(II).C(:, IR) = DD.cR(:, DD.lsx(h, 1));
            end
            II = i + TT;
            CE(II) = CE(II) + 1;
            MSS(II).F(CE(II)) = DD.O;
            MSS(II).E(CE(II)) = DD.eO;
            MSS(II).C(:, CE(II)) = DD.cO;
            MSS(II).L(:, CE(II)) = DD.LO;
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, S, VAL] = CALCULS_WORKER(Const, i, MSS, VAL, PH)
%--------------------------------------------------------------------------
% Create and allocating empty sets
A.cO = MSS(i).C(:, PH);
A.LO = MSS(i).L(:, PH);
A.eO = MSS(i).E(PH);
A.O = MSS(i).F(PH);
max_L = max(A.LO);
DIVIS = find(A.LO == max_L);
DELTA  = 2*max_L/3;
S = length(DIVIS);
[A.L, A.R, A.eL, A.eR] = deal(zeros(1, S));
A.lL = A.LO*ones(1, S);
A.cL = A.cO*ones(1, S); A.cR = A.cL;

for g = 1:S
    A.cL(DIVIS(g), g) = A.cL(DIVIS(g), g) - DELTA;
    om_point = abs(VAL.b - VAL.a).*A.cL(:,g) + VAL.a;
    A.L(g) = feval(Const.f, om_point);
    A.cR(DIVIS(g), g) = A.cR(DIVIS(g), g) + DELTA;
    om_point = abs(VAL.b - VAL.a).*A.cR(:,g) + VAL.a;
    A.R(g) = feval(Const.f, om_point);
end
[~, A.lsx] = sort([min(A.L, A.R)' DIVIS], 1);
for g = 1:S
    A.lL(DIVIS(A.lsx(1:g, 1)), g) = DELTA/2;
    A.eL(g) = VAL.fcount + 1;
    A.eR(g) = VAL.fcount + 2;
    VAL.fcount = A.eR(g);
end
A.LO = A.lL(:, S);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
%--------------------------------------------------------------------------
function [VAL, SS] = Arewedone(minval, VAL, SS)
%--------------------------------------------------------------------------
if (abs(minval - VAL.fMinBeforeImpr) > VAL.Yaro_epsilon*abs(minval))
    VAL.fMinNotImpr     = 0;
    VAL.fMinBeforeImpr  = minval;
end

% Show iteration stats
if SS.showITS == 1
    VAL.time = toc;
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
hulls    = cell(3, index);
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
if (abs(minval - MV(i_min)) < 0.01)
    if (abs(minval - VAL.fMinBeforeImpr) < VAL.Yaro_epsilon*abs(minval))
        VAL.fMinNotImpr     = VAL.fMinNotImpr + 1;
    else
        VAL.fMinNotImpr     = 0;
        VAL.fMinBeforeImpr  = minval;
    end
end
if (VAL.fMinNotImpr > 5) && (mod(VAL.fMinNotImpr, 5) ~= 0)
    i_m = i_m - round(6*(i_m)/8);
end

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
            TT = find(MSS(i).F(hulls{1, i}) <= slope*d(i) + const + 1E-12);
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