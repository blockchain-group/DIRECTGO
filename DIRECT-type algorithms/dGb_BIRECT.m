function [minima, xatmin, history] = dGb_BIRECT(Problem, bounds, opts)
%--------------------------------------------------------------------------
% Function   : dGb_BIRECT
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dGb_BIRECT(Problem, bounds, opts)
%       d         - dynamic memory management in data structure
%       BIRECT    - BIRECT(BIsecting RECTangles) algorithm
%       Gb        - globally-biased
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
% Paulavicius, R., Chiter, L., and Zilinskas, J. "Global optimization 
% based on bisection of rectangles, function values at diagonals, and a
% set of Lipschitz constants. Journal of Global Optimization. (2018). 
% DOI 10.1007/s10898-016-0485-6 
% 
% Paulavicius, R., Sergeyev, Y. D., Kvasov, D. E., & Zilinskas, J. 
% "Globally-biased BIRECT algorithm with local accelerators for expensive 
% global optimization" Expert Systems with Applications. (2020).
% DOI 10.1016/j.eswa.2019.113052
%--------------------------------------------------------------------------

% Get options
SS = Options(opts, nargout);

% Alocate sets and create initial variables
[MSS, CE, VAL] = Alocate(bounds, SS);

% Initialization step
[MV, MSS, CE, VAL, minval, xatmin] = Initialization(Problem,...
    MSS, CE, VAL, SS);

while VAL.perror > SS.TOL                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    POH = Find_po(MSS, CE, SS, minval, MV(1, :), VAL);

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
    history   = VAL.history(1:(VAL.itctr - 1), 1:4);
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
function SS = Options(opts,narg)
%--------------------------------------------------------------------------
% Determine option values
if nargin < 3 && (isempty(opts))
    opts = [];
end
getopts(opts, 'maxits', 100, 'maxevals', 100000, 'maxdeep', 1000,...
    'testflag', 0, 'globalmin', 0, 'tol', 0.01, 'showits', 1,...
    'ep', 1e-4);

SS.G_nargout = narg;     % output arguments
SS.MAXits    = maxits;   % maximum of iterations
SS.MAXevals  = maxevals; % maximum # of function evaluations
SS.MAXdeep   = maxdeep;  % maximum number of side divisions
SS.showITS   = showits;  % print iteration stat
SS.TOL       = tol;      % allowable relative error if f_reach is set
SS.TESTflag  = testflag; % terminate if within a relative tolerence of f_opt
SS.globalMIN = globalmin;% minimum value of function
SS.ep        = ep;       % global/local weight parameter
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = Alocate(bounds, SS)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
VAL.a          = bounds(:, 1);              % left bound
VAL.b          = bounds(:, 2);              % right bound
VAL.n          = size(bounds, 1);           % dimension
VAL.time       = 0;                         % initial time
VAL.fcount     = 2;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
CE             = zeros(1, VAL.n*SS.MAXdeep);% collum counter

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'L', zeros(VAL.n, 1), 'RECT', []);

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
MSS(1).RECT(1).p1 = ones(VAL.n, 1)./3;             
MSS(1).RECT(1).p2 = ones(VAL.n, 1).*(2/3);               
MSS(1).RECT(1).c = (MSS(1).RECT(1).p1 + MSS(1).RECT(1).p2)/2;         
MSS(1).RECT(1).f1 = feval(Problem.f, abs(VAL.b -...
    VAL.a).*MSS(1).RECT(1).p1 + VAL.a);
MSS(1).RECT(1).f2 = feval(Problem.f, abs(VAL.b -...
    VAL.a).*MSS(1).RECT(1).p2 + VAL.a);

if MSS(1).RECT(1).f1 < MSS(1).RECT(1).f2
    MSS(1).RECT(1).fmin = MSS(1).RECT(1).f1;
    MSS(1).RECT(1).pmin = MSS(1).RECT(1).p1;
else
    MSS(1).RECT(1).fmin = MSS(1).RECT(1).f2;
    MSS(1).RECT(1).pmin = MSS(1).RECT(1).p2;
end
[Fmin, MV(1), MSS(1).F(1), VAL.fMinBeforeImpr] = deal(MSS(1).RECT(1).fmin);
Xmin = MSS(1).RECT(1).pmin;

MSS(1).L(:, 1) = ones(VAL.n, 1);
MSS(1).Diam = roundn((2/3)*sqrt(sum(MSS.L(:, 1).^2)), -12);
VAL.fMinNotImpr = 0;

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
                [MSS(i).L, MSS(i).F, MSS(i).RECT] = deal([]);
            end
        else
            C = 1:CE(i); C(POH{i}) = [];
            pp = min(POH{i}):length(C);
            MSS(i).L(:, pp) = MSS(i).L(:, C(pp));
            MSS(i).F(pp) = MSS(i).F(C(pp));
            MSS(i).RECT(:, pp) = MSS(i).RECT(:, C(pp));
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
            I = i + 1;
            if I > length(CE), CE(I) = 0; end
            if CE(I) == 0
                MSS(I).L = zeros(VAL.n, 100);
                MSS(I).F = zeros(1, 100);
            end
            % Initial calculations       
            max_L = max(MSS(i).L(:, PH{i}(j)));
            lsa = find(MSS(i).L(:, PH{i}(j)) == max_L);
            ls = lsa(1);
            dl = max_L/2;
            IL = CE(I) + 1;             
            IR = CE(I) + 2;             
            CE(I) = IR;
            [MSS(I).L(:, IL), MSS(I).L(:, IR)] = deal(MSS(i).L(:, PH{i}(j)));
            [MSS(I).L(ls, IL), MSS(I).L(ls, IR)] = deal(dl);
            [MSS(I).RECT(IL), MSS(I).RECT(IR)] = deal(MSS(i).RECT(PH{i}(j)));
    
            % New sampling points
            if MSS(I).RECT(IL).p1(ls) < MSS(I).RECT(IL).p2(ls)
                MSS(I).RECT(IR).p1(ls) = MSS(I).RECT(IL).p1(ls) + dl;
                MSS(I).RECT(IL).p2(ls) = MSS(I).RECT(IR).p2(ls) - dl;
            else
                MSS(I).RECT(IR).p1(ls) = MSS(I).RECT(IL).p1(ls) - dl;
                MSS(I).RECT(IL).p2(ls) = MSS(I).RECT(IR).p2(ls) + dl;
            end
            
            % Update center points for both rectangles
            MSS(I).RECT(IL).c = (MSS(I).RECT(IL).p1 + MSS(I).RECT(IL).p2)/2;
            MSS(I).RECT(IR).c = (MSS(I).RECT(IR).p1 + MSS(I).RECT(IR).p2)/2;
            
            % Evaluate function at points
            MSS(I).RECT(IR).f1 = feval(Problem.f, abs(VAL.b -...
                VAL.a).*MSS(I).RECT(IR).p1 + VAL.a);
            MSS(I).RECT(IL).f2 = feval(Problem.f, abs(VAL.b -...
                VAL.a).*MSS(I).RECT(IL).p2 + VAL.a);
            
            % Find/update the minimal obj. value and the minimum point
            if MSS(I).RECT(IL).f1 < MSS(I).RECT(IL).f2
                MSS(I).RECT(IL).fmin = MSS(I).RECT(IL).f1;
                MSS(I).RECT(IL).pmin = MSS(I).RECT(IL).p1;
            else
                MSS(I).RECT(IL).fmin = MSS(I).RECT(IL).f2;
                MSS(I).RECT(IL).pmin = MSS(I).RECT(IL).p2;
            end
            
            if MSS(I).RECT(IR).f1 < MSS(I).RECT(IR).f2
                MSS(I).RECT(IR).fmin = MSS(I).RECT(IR).f1;
                MSS(I).RECT(IR).pmin = MSS(I).RECT(IR).p1;
            else
                MSS(I).RECT(IR).fmin = MSS(I).RECT(IR).f2;
                MSS(I).RECT(IR).pmin = MSS(I).RECT(IR).p2;
            end
            
            % Update min func. value for POH(1) rectangle
            MSS(I).F(IL) = MSS(I).RECT(IL).fmin;
            % Add new min value of the new rectangle
            MSS(I).F(IR) = MSS(I).RECT(IR).fmin;
            
            MSS(I).Diam = roundn((2/3)*sqrt(sum(MSS(I).L(:, IL).^2)), -12);
            VAL.fcount = VAL.fcount + 2;
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
if (abs(minval - VAL.fMinBeforeImpr) < 0.01*abs(minval))
    VAL.fMinNotImpr = VAL.fMinNotImpr + 1;
else
    VAL.fMinNotImpr    = 0;
    VAL.fMinBeforeImpr = minval;
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
    VAL.history(VAL.itctr,1) = VAL.itctr;
    VAL.history(VAL.itctr,2) = VAL.fcount;
    VAL.history(VAL.itctr,3) = minval;
    VAL.history(VAL.itctr,4) = VAL.time;
end
% Update iteration number
VAL.itctr = VAL.itctr + 1;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function S = Find_po(MSS, CE, SS, minval, MV, VAL)
% Find all rects on hub
index    = size(MSS, 2);
d        = [MSS.Diam];
hulls    = cell(3, index);
hs       = nan(2, index);

for i = 1:index
    if CE(i) ~= 0
        ties     = find(MSS(i).F(1:CE(i)) - MV(i) <= 1E-12);
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
i_m   = i_min - jj;
if (VAL.fMinNotImpr > 2*VAL.n) && mod(VAL.fMinNotImpr, 10) ~= 0
    i_m = i_m - floor(6*(i_m)/8);
end

jl        = i_m + jj;
d_min = MV(1, :);
S_1   = cell(1, index);
S_1(1:jl)   = hulls(1, 1:jl);
Su    = cell(4, index);

if i_m > 1
    a1 = d(i_min);
    b1 = MSS(i_min).F(hs(2, i_min));
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
xatmin = MSS(Least).RECT(MV(2, Least)).pmin;
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