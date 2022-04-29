function [minima, xatmin, history] = dADC(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : dADC
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dADC(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
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
[SS, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[MSS, CE, VAL, third] = Alocate(SS, VAL);

% Initialization step
[MV, MSS, CE, VAL, minval, xatmin] = Initialization(Problem,...
    MSS, CE, VAL, SS);

while VAL.perror > SS.TOL                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    [POH, VAL] = Find_poh(MSS, CE, SS, minval, MV(1, :), VAL);
    
    % Subdivide potential optimalhyper-rectangles
    [VAL, MSS, CE, MV, minval, xatmin] = Calulcs(VAL, Problem,...
        MSS, CE, POH, minval, xatmin, third);

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
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL, third] = Alocate(SS, VAL)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
[VAL.time, VAL.global] = deal(0);           % initial time
VAL.fcount     = 1;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
CE             = zeros(1, VAL.n*SS.MAXdeep);% collum counter

third       = zeros(1, SS.MAXdeep);         % delta values
third(1)    = 1/3;                          % first delta
for i = 2:SS.MAXdeep                        % all delta
    third(i)  = (1/3)*third(i - 1);
end

% alociate MAIN sets
MSS = struct('F', zeros(1), 'CC', zeros(VAL.n, SS.MAXevals),...
    'FF', zeros(1), 'L', zeros(VAL.n, 1), 'RECT', []);

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
CE(1) = 1;
MSS(1).RECT(1).p1 = zeros(VAL.n, 1);             
MSS(1).RECT(1).p2 = ones(VAL.n, 1);                     
MSS(1).RECT(1).f1 = feval(Problem.f, VAL.a);
MSS(1).RECT(1).f2 = feval(Problem.f, VAL.b);
MSS(1).L(:, 1) = zeros(VAL.n, 1);
MSS(1).CC(:, 1) = zeros(VAL.n, 1);
MSS(1).CC(:, 2) = ones(VAL.n, 1);
VAL.e = 2;
Fmin = min(MSS(1).RECT(1).f1, MSS(1).RECT(1).f2);

if MSS(1).RECT(1).f1 < MSS(1).RECT(1).f2
    Xmin = MSS(1).RECT(1).p1;
else
    Xmin = MSS(1).RECT(1).p2;
end
[MSS(1).FF(1), MV(1)] = deal((MSS(1).RECT(1).f1 + MSS(1).RECT(1).f2)/2);
[MSS(1).F(1), VAL.fMinBeforeImpr] = deal(Fmin);
VAL.fMinNotImpr = 0;
MSS(1).Diam = 1/2*norm((1/3*(ones(VAL.n, 1))).^(MSS(1).L(:, 1)), 2);

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
    VAL.history(VAL.itctr, 2) = VAL.e;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [VAL, MSS, CE, MV, Fmin, Xmin] = Calulcs(VAL,...
    Problem, MSS, CE, POH, Fmin, Xmin, third)
%--------------------------------------------------------------------------
[MSS, CE, VAL, Fmin, Xmin] = STORAGE(Problem, MSS, CE, VAL, POH,...
    Fmin, Xmin, third);

for i = size(POH, 2):-1:1                
    if ~isempty(POH{i})         
        if (CE(i) - size(POH{i}, 2)) == 0
            if find(CE ~= 0, 1, 'first') == i
                [MSS(i).L, MSS(i).FF, MSS(i).F, MSS(i).RECT,...
                    MSS(i).diam] = deal([]);
            end
        else
            C = 1:CE(i); C(POH{i}) = [];
            pp = min(POH{i}):length(C);
            MSS(i).L(:, pp) = MSS(i).L(:, C(pp));
            MSS(i).F(pp) = MSS(i).F(C(pp));
            MSS(i).FF(pp) = MSS(i).FF(C(pp));
            MSS(i).RECT(:, pp) = MSS(i).RECT(:, C(pp));
        end
        CE(i) = CE(i) - size(POH{i}, 2);
    end
end

% Find minima values
MV = Find_min(CE, MSS);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [M, CE, VAL, Fmin, Xmin] = STORAGE(Problem, M, CE, VAL,...
    PH, Fmin, Xmin, third)
%--------------------------------------------------------------------------
for i = 1:size(M, 2)
    if ~isempty(PH{i})
        for j = 1:size(PH{i}, 2)
            I = i + 1;
            if I > size(CE, 2)
                CE(I) = 0;    
            end
            if CE(I) == 0
                M(I).L = zeros(VAL.n, 100);
                M(I).F = zeros(1, 100);
                M(I).FF = zeros(1, 100);
            end
            IL = CE(I) + 1; IR = CE(I) + 2; IO = CE(I) + 3; CE(I) = IO; 
            O = PH{i}(j);
            % Initial calculations
            ls = find(M(i).L(:, O) == min(M(i).L(:, O)), 1, 'first');
            delta = third(min(M(i).L(:, O)) + 1);
            M(i).L(ls, O) = M(i).L(ls, O) + 1;
            [M(I).L(:, IL), M(I).L(:, IR), M(I).L(:, IO)] = deal(M(i).L(:, O));
            [M(I).RECT(IL), M(I).RECT(IR)] = deal(M(i).RECT(O));
            
            % New sampling points
            if M(i).RECT(O).p1(ls) < M(i).RECT(O).p2(ls)
                M(I).RECT(IL).p1(ls) = M(I).RECT(IL).p1(ls) + 2*delta;
                M(I).RECT(IL).p2(ls) = M(I).RECT(IL).p2(ls) - 2*delta;
                M(I).RECT(IR).p1 = M(I).RECT(IL).p1; 
                M(i).RECT(O).p2 = M(I).RECT(IL).p2;
            else
                M(I).RECT(IL).p1(ls) = M(I).RECT(IL).p1(ls) - 2*delta;
                M(I).RECT(IL).p2(ls) = M(I).RECT(IL).p2(ls) + 2*delta;
                M(I).RECT(IR).p1 = M(I).RECT(IL).p1; 
                M(i).RECT(O).p2 = M(I).RECT(IL).p2;
            end
            
            % Update center points for both rectangles
            [M(I).RECT(IR).f1, M(I).RECT(IL).f1] = deal(feval(Problem.f,...
                (abs(VAL.b - VAL.a).*M(I).RECT(IL).p1 + VAL.a)));
            [M(I).RECT(IL).f2, M(i).RECT(O).f2] = deal(feval(Problem.f,...
                (abs(VAL.b - VAL.a).*M(I).RECT(IL).p2 + VAL.a)));
            
            M(I).FF(IO) = (M(i).RECT(O).f1 + M(i).RECT(O).f2)/2; 
            M(I).F(IO) = min(M(i).RECT(O).f1, M(i).RECT(O).f2);      
            M(I).FF(IL) = (M(I).RECT(IL).f1 + M(I).RECT(IL).f2)/2; 
            M(I).F(IL) = min(M(I).RECT(IL).f1, M(I).RECT(IL).f2);
            M(I).FF(IR) = (M(I).RECT(IR).f1 + M(I).RECT(IR).f2)/2; 
            M(I).F(IR) = min(M(I).RECT(IR).f1, M(I).RECT(IR).f2);

            idx = find(~any(bsxfun(@minus, M(1).CC(:, 1:VAL.e), M(I).RECT(IL).p1)), 1, 'first');
            
            if isempty(idx)
                VAL.e = VAL.e + 1;
                M(1).CC(:, VAL.e) = M(I).RECT(IL).p1;
            end
            idx = find(~any(bsxfun(@minus, M(1).CC(:, 1:VAL.e), M(I).RECT(IL).p2)), 1, 'first');
            if isempty(idx)
                VAL.e = VAL.e + 1;
                M(1).CC(:, VAL.e) = M(I).RECT(IL).p2;
            end

            f_calls = [M(i).RECT(O).f1, M(i).RECT(O).f2,...
                M(I).RECT(IL).f1, M(I).RECT(IL).f2,...
                M(I).RECT(IR).f1, M(I).RECT(IR).f2];
            pts = [M(i).RECT(O).p1, M(i).RECT(O).p2,...
                M(I).RECT(IL).p1, M(I).RECT(IL).p2,...
                M(I).RECT(IR).p1, M(I).RECT(IR).p2];
            [ff, index] = min(f_calls);
    
            if ff < Fmin, Fmin = ff; Xmin = pts(:, index); end

            M(I).Diam = round(1/2*norm((1/3*(ones(VAL.n,...
                1))).^(M(I).L(:, IL)), 2), 16);
            VAL.fcount = VAL.fcount + 2;
            M(I).RECT(IO) = M(i).RECT(O);
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
if (abs(minval - VAL.fMinBeforeImpr) < 0.01*abs(minval))
    VAL.fMinNotImpr = VAL.fMinNotImpr + 1;
else
    VAL.fMinNotImpr    = 0;
    VAL.fMinBeforeImpr = minval;
end

VAL.time = toc;

if SS.showITS == 1
    fprintf(...
    'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, minval, VAL.time, VAL.e);
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
if VAL.e > SS.MAXevals
    disp('Exceeded max fcn evals. Increase maxevals'); VAL.perror = -1;
end

% Store History
if SS.G_nargout == 3
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.e;
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
function [S, VAL] = Find_poh(MSS, CE, SS, minval, MV, VAL)
% Find all rects on hub
index    = size(MSS, 2);
d        = [MSS.Diam];
hulls    = cell(3, index);
hs       = nan(2, index);

for i = 1:index
    if CE(i) ~= 0
        ties     = find(MSS(i).FF(1:CE(i)) - MV(i) <= 1E-12);
        if ~isempty(ties)
            hulls{1, i} = ties;
            hulls{2, i} = d(i)*ones(1, size(ties, 2));
            hulls{3, i} = MSS(i).FF(ties);
            hulls{4, i} = i*ones(1, size(ties, 2));
            [hs(1, i), hs(2, i)] = min((MSS(i).FF(1:CE(i)) - minval +...
                max(SS.ep*abs(minval), 1E-8))/d(i));
        end
    end
end

jj    = find(~cellfun(@isempty, hulls(1, :)), 1, 'first');
[i_min, jl] = deal(find(hs(1, :) == min(hs(1, :))));

if VAL.fMinNotImpr > VAL.n + 1  
    [~, findx] = min(MV(1, :));
    q   = sum(MSS(find(CE ~= 0, 1, 'first')).L(:, 1));
    p   = sum(MSS(findx).L(:, 1));
    jl = fix((q + p)/2) + 1;
    VAL.global = VAL.global + 1;
    if VAL.global > 2^(VAL.n + 1)
        [VAL.global, VAL.fMinNotImpr] = deal(0);
    end
end

d_min = MV(1, :);
S_1   = cell(1, index);
S_1(1:jl)   = hulls(1, 1:jl);
Su    = cell(4, index);

if jl - jj > 1
    a1 = d(i_min);
    b1 = MSS(i_min).FF(hs(2, i_min));
    a2 = d(jj);
    b2 = d_min(jj);
    % The line is defined by: y = slope*x + const
    slope = (b2 - b1)/(a2 - a1);
    const = b1 - slope*a1;
    for i = 1:jl
        if ~isempty(hulls{3, i})
            TT = find(hulls{3, i} <= slope*d(i) + const + 1E-8);
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
        S_2 = cell2mat(hulls(1:4, 1:jl));
        S_2 = sortrows(S_2.', 2).';
    else
        S_2 = sortrows(S_2.', 2).';
    end
    
    % S_2 now contains all points in S_1 which lies on or below the line
    % Find the points on the convex hull defined by the points in S_2
    xx = S_2(2, :);
    yy = S_2(3, :);
    h = conhull(xx, yy); % conhull is an internal subfunction
    invers = S_2(:, h);
    S = cell(1, index);
    for i = 1:jl
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
function MV = Find_min(CE, MSS)
%--------------------------------------------------------------------------
MV = nan(2, size(MSS, 2));

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        [MV(1, i), MV(2, i)] = min(MSS(i).FF(1:CE(i)));
    end
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