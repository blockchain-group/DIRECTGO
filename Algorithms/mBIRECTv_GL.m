function [minima, xatmin, history] = mBIRECTv_GL(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : mBIRECTv_GL
% Written by : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Written by : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 05/30/2023
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = mBIRECTv_GL(Problem, opts, bounds)
%       m       - Mapping Technique for handling linear constraints
%       BIRECTv - Hyper-rectangular partitioning based on 1-Dimensional 
%                  Bisection and objective function evaluations at two
%                  points
%       G       - Enhancing the global search
%       L       - Enhancing the local search
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f          = Objective function handle
%                 Problem.constraint = Constraint function handle
%
%       opts    - MATLAB structure which contains options.
%                 opts.maxevals  = max. number of function evals
%                 opts.maxits    = max. number of iterations
%                 opts.globalmin = globalmin (if known)
%                 opts.globalxmin = globalxmin (if known)
%                 opts.testflag  = 1 if globalmin known, 0 otherwise
%                 opts.dimension = problem dimension
%                 opts.showits   = 1 print iteration status
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
% Partitioning strategy taken from:
%--------------------------------------------------------------------------
% Chiter, L. Experimental Data for the preprint "Diagonal Partitioning 
% Strategy Using Bisection of Rectangles and a Novel Sampling 538
% Scheme", 2023. https://doi.org/10.17632/x9fpc9w7wh.2.
%
% Selection of potential optimal hyper-rectangles taken from:
%--------------------------------------------------------------------------
% Stripinis, L., Paulavicius, R., Zilinskas, J.: Improved scheme for
% selection of potentially optimal hyperrectangles in DIRECT. Optimization
% Letters (2018). ISSN 1862-4472, 12 (7), 1699-1712,
% DOI: 10.1007/s11590-017-1228-4
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end

% Get options
[OPTI, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[VAL, MSS, RECT] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, Xmin, Fmin, MSS, RECT] = Initialization(VAL, OPTI,...
    Problem, MSS, RECT);

while VAL.perror > OPTI.TOL                                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    POH = Selection(VAL, MSS, Xmin);

    % Subdivide potential optimalhyper-rectangles
    [MSS, VAL, RECT] = Subdivision(VAL, Problem, MSS, POH, RECT);
    
    % Update minima and check stopping conditions
    [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS, RECT);
end                                                         % End of while

% Return value
minima        = Fmin;
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
% Function  : Selection
% Purpose   : Selection of potential optimal hyper-rectangles
%--------------------------------------------------------------------------
function POH = Selection(VAL, MSS, Xmin)
%--------------------------------------------------------------------------
% Calculate Euclidean Distatnces
Euclid_dist = sum((Xmin(:, 1) - MSS.CC(:, 1:VAL.I)).^2, 1).^0.5;

% Identify potential optimal hyper-rectangles
S = Find_poh(MSS.FF(1:VAL.I), MSS.DD(1:VAL.I), 1:VAL.I);

D = Find_poh(Euclid_dist(1:VAL.I), MSS.DD(1:VAL.I), 1:VAL.I);

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
function boxes = Find_poh(fc, szes, indexs)
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
getopts(opts, 'maxits', 1000,'maxevals', 100000, 'testflag', 0, ...
    'tol', 0.01, 'showits', 1, 'dimension', 1, 'globalmin', 0, 'globalxmin', 0);

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

VAL = FindVertices([VAL.a, VAL.b], Problem, VAL);
[A, B] = FindVerticesolds([VAL.a, VAL.b], Problem);
V = con2vert(A, B);

for i = 1:VAL.n
    VAL.a(i) = max([min(V(:, i)), VAL.a(i)]);
    VAL.b(i) = min([max(V(:, i)), VAL.b(i)]);
    V((V(:, i) <= VAL.a(i)), i) = VAL.a(i);
    V((V(:, i) >= VAL.b(i)), i) = VAL.b(i);
end
VAL.C = mean(V, 1)';
OPTI.G_nargout = narg;     % output arguments
OPTI.MAXits    = maxits;   % Fmax of iterations
OPTI.MAXevals  = maxevals; % Fmax # of function evaluations
OPTI.TESTflag  = testflag; % terminate if global minima is known
OPTI.showITS   = showits;  % print iteration stat
OPTI.TOL       = tol;      % allowable relative error if f_reach is set
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
VAL.perror  = 10;                       % initial perror

% alociate MAIN sets
z   = round(OPTI.MAXevals);
MSS = struct('FF', zeros(1, (z/2)),'FFF', zeros(1, (z/2)), 'CCC', -ones(VAL.n, (z/2)),...
    'CC', -ones(VAL.n, (z/2)), 'DD', -ones(1, (z/2)), 'LL', zeros(VAL.n, (z/2)));

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
[VAL.I, VAL.e, VAL.itctr] = deal(1);      
RECT(1).p1 = ones(VAL.n, 1).*(1/3);       
RECT(1).p2 = ones(VAL.n, 1);              

RECT(1).f1 = feval(Problem.f, Topology(abs(VAL.b - VAL.a).*RECT(1).p1 + VAL.a, VAL));
RECT(1).f2 = feval(Problem.f, Topology(abs(VAL.b - VAL.a).*RECT(1).p2 + VAL.a, VAL));

if RECT(1).f1 < RECT(1).f2
    [Fmin, MSS.FF(1)] = deal(RECT(1).f1); Xmin = RECT(1).p1;
else
    [Fmin, MSS.FF(1)] = deal(RECT(1).f2); Xmin = RECT(1).p2;
end

MSS.FFF(1) = RECT(1).f2;
MSS.CCC(:, 1) = RECT(1).p2;

MSS.CC(:, 1) = ones(VAL.n, 1)./2;
MSS.LL(:, 1) = ones(VAL.n, 1);
MSS.DD(1) = (2/3)*sqrt(sum(MSS.LL(:, 1).^2));

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
    VAL.history(VAL.itctr, 2) = VAL.I + VAL.e;
    VAL.history(VAL.itctr, 3) = Fmin;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
%--------------------------------------------------------------------------
function [VAL, Fmin, Xmin] = Arewedone(OPTI, VAL, MSS, RECT)
%--------------------------------------------------------------------------
[Fmin, fminindex] =  min(MSS.FF(1:VAL.I));
if RECT(fminindex).f1 < RECT(fminindex).f2
    Xmin = RECT(fminindex).p1;
else
    Xmin = RECT(fminindex).p2;
end
VAL.time = toc;

if OPTI.showITS == 1                % Show iteration stats
    fprintf(...
    'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, Fmin, VAL.time, VAL.I + VAL.e);
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

if VAL.I + VAL.e > OPTI.MAXevals       % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.I + VAL.e;
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
function [MSS, VAL, RECT] = Subdivision(VAL, O, MSS, POH, RECT)
%--------------------------------------------------------------------------
while isempty(POH) == 0
    lsa = find(MSS.LL(:, POH(1)) ==  max(MSS.LL(:, POH(1))));
    for h = 1:length(lsa)
        VAL.I = VAL.I + 1; 
        delta = (max(MSS.LL(:, POH(1))))/2;
        ls = find(MSS.LL(:, POH(1)) == max(MSS.LL(:, POH(1))), 1, "first");
        MSS.LL(:, VAL.I) = MSS.LL(:, POH(1)); 
        RECT(VAL.I) = RECT(POH(1));
        [MSS.LL(ls, VAL.I), MSS.LL(ls, POH(1))] = deal(delta);
        MSS.CC(:, VAL.I) = MSS.CC(:, POH(1));
        if RECT(POH(1)).p1(ls) < RECT(POH(1)).p2(ls)
            MSS.CC(ls, POH(1)) = MSS.CC(ls, POH(1)) - delta/2;
            MSS.CC(ls, VAL.I) =  MSS.CC(ls, VAL.I) + delta/2;
            RECT(VAL.I).p1(ls) = RECT(POH(1)).p1(ls) + delta*(2/3);
            RECT(POH(1)).p2(ls) = RECT(VAL.I).p2(ls) - delta*2;
            idx = find(~any(bsxfun(@minus, MSS.CCC(:, 1:VAL.e), RECT(POH(1)).p2)), 1, 'first');
            if isempty(idx)
                VAL.e = VAL.e + 1;
      	        RECT(POH(1)).f2 = feval(O.f, Topology(abs(VAL.b - VAL.a).*RECT(POH(1)).p2 + VAL.a, VAL));
                MSS.CCC(:, VAL.e) = RECT(POH(1)).p2;
                MSS.FFF(VAL.e) = RECT(POH(1)).f2;
            else
                RECT(POH(1)).f2 = MSS.FFF(idx);
            end
        else
            MSS.CC(ls, POH(1)) = MSS.CC(ls, POH(1)) + delta/2;
            MSS.CC(ls, VAL.I) =  MSS.CC(ls, VAL.I) - delta/2;
            RECT(VAL.I).p1(ls) = RECT(POH(1)).p1(ls) - delta*(2/3);
            RECT(POH(1)).p2(ls) = RECT(VAL.I).p2(ls) + delta*2;
            idx = find(~any(bsxfun(@minus, MSS.CCC(:, 1:VAL.e), RECT(POH(1)).p2)), 1, 'first');
            if isempty(idx)
                VAL.e = VAL.e + 1;
                RECT(POH(1)).f2 = feval(O.f, Topology(abs(VAL.b - VAL.a).*RECT(POH(1)).p2 + VAL.a, VAL));
                MSS.CCC(:, VAL.e) = RECT(POH(1)).p2;
                MSS.FFF(VAL.e) = RECT(POH(1)).f2;
            else
                RECT(POH(1)).f2 = MSS.FFF(idx);
            end
        end
        RECT(VAL.I).f1 = feval(O.f, Topology(abs(VAL.b - VAL.a).*RECT(VAL.I).p1 + VAL.a, VAL));
        MSS.FF(POH(1)) = min([RECT(POH(1)).f1, RECT(POH(1)).f2]);
        MSS.FF(VAL.I) = min([RECT(VAL.I).f1, RECT(VAL.I).f2]);
        [MSS.DD(VAL.I), MSS.DD(POH(1))] = deal((2/3)*sqrt(sum(MSS.LL(:, POH(1)).^2)));
        if MSS.FF(POH(1)) > MSS.FF(VAL.I)
            POH(1) = VAL.I;
        end
    end
    POH(1) = [];
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  FindVertices
% Purpose    :  Cover hyper-rectangle by simplex
%--------------------------------------------------------------------------
function VAL = FindVertices(bounds, Problem, VAL)
% Find intersecting vertices accorting to linear constrains
warning('off','all');
B = -feval(Problem.constraint, zeros(1, size(bounds, 1)))';
Idt = repelem(eye(size(bounds, 1)), 1, 1);
A = zeros(length(B), VAL.n);
for i = 1:VAL.n
    A(:, i) = feval(Problem.constraint, Idt(:, i))' + B;
end
VAL.AA = A;
VAL.BB = B;
return

%--------------------------------------------------------------------------
% Function   :  FindVertices
% Purpose    :  Cover hyper-rectangle by simplex
%--------------------------------------------------------------------------
function [A, B] = FindVerticesolds(bounds,Problem)
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
return

%--------------------------------------------------------------------------
% Function   :  con2vert
% Purpose    :  convert a convex set of constraint inequalities into the set
%               of vertices at the intersections of those inequalities
% Source     :  Stephen Becker (2022). CON2VERT - constraints to vertices, 
%               updated (https://www.mathworks.com/matlabcentral/fileexchange
%               /91595-con2vert-constraints-to-vertices-updated), MATLAB 
%               Central File Exchange. Retrieved September 9, 2022.
%--------------------------------------------------------------------------
function V = con2vert(A, b)
c = A\b; 
if ~all(A*c < b)
    dim = size(A,2);
    m   = size(A,1);
    opts = optimoptions('linprog','ConstraintTolerance', 1e-8);
    [x_and_t] = linprog([zeros(dim,1); 1], [A,-ones(m, 1)], b, [], [], [], [], opts);
    c = x_and_t(1:dim);
end
b = b - A*c;    
D = bsxfun( @times, A, 1./b );
k = convhulln(D);
G = zeros(size(k, 1),size(D,2));
o = ones( size(k, 2), 1 );  

badIndices = [];
for ix = 1:size(k, 1)
    F = D(k(ix, :), :);
    G(ix, :)=F\o;
    [~, LASTID] = lastwarn;
    if any( isnan(G(ix,:)) ) || strcmpi(LASTID,'MATLAB:singularMatrix')
        badIndices(end + 1) = ix; %#ok<AGROW> 
    end
    lastwarn('');  
end
G(badIndices, :) = [];
V1 = round( G*randn(size(G,2),1),  10, 'significant' );
V2 = round( G*randn(size(G,2),1),  10, 'significant' );
[~,I1]=unique(V1,'stable');
[~,I2]=unique(V2,'stable');
I = union(I1, I2);
G = G(I, :);
V = round(bsxfun(@plus, G, c'), 8);      
return

%--------------------------------------------------------------------------
% Function   :  Find topology point
% Purpose    :  Find alternative point in the domain
%--------------------------------------------------------------------------
function point_t = Topology(point, VAL)
%--------------------------------------------------------------------------
point_t = point;
% Find direction of the point
direction = VAL.C - point;
direction(direction < 0) = -1;
direction(direction >= 0) = 1;
plane = zeros(VAL.n, 1);
for i = 1:VAL.n
    if direction(i) == -1
        plane(i) = VAL.b(i);
    else
        plane(i) = VAL.a(i);
    end
end

% Find a parameterization of the line and intersaction on boundary
betha = point - VAL.C;
betha(betha == 0) = 1e-16; % Add small error
t = (plane - point)./betha;
point_b = round(point + betha.*t, 10);
bound_index = find(point_b == plane);
for i = 1:length(bound_index)
    point_b = round(point + betha*t(bound_index(i)), 10);
    if all(point_b <= VAL.b) && all(point_b >= VAL.a)
        break
    end
end

% Check feasibility of boundary point
constraints = sum(VAL.AA.*point_b', 2) - VAL.BB;
if_i = find(constraints > 0);
if ~isempty(if_i)
    point_t = zeros(VAL.n, length(if_i));
    for i = 1:length(if_i)
        t_t = sum(betha.*VAL.AA(if_i(i), :)');
        value = VAL.BB(if_i(i), :) - sum(point.*VAL.AA(if_i(i), :)');
        t = value/t_t;
        point_t(:, i) = point + betha.*t;
    end
    E_dist = sum((VAL.C - point_t).^2, 1).^0.5;
    if isnan(E_dist)
        E_dist = 0;
    end
    ee = find(E_dist == min(E_dist), 1);
    point_t = point_t(:, ee);
    E_b = sum((VAL.C - point_b).^2, 1).^0.5;
    E_d = sum((VAL.C - point_t).^2, 1).^0.5;
    r_ratio = (E_d./E_b);
    point_t = VAL.C + betha.*r_ratio;
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