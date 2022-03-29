function [minima, xatmin, history] = dDirect_GLce(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : dDirect_GLce
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for problems with various
%              constraints constraints
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dDirect_GLce(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       G         - Enhancing the global search
%       L         - Enhancing the local search
%       c         - Constraint handling
%       e         - Additional tolerance for constraint handling
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f           = Objective function handle
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
%                 opts.ept       = tollerance for constrains
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
% Selection of potential optimal hyper-rectangles taken from:
%--------------------------------------------------------------------------
% Stripinis, L., Paulavicius, R., Zilinskas, J.: Improved scheme for
% selection of potentially optimal hyperrectangles in DIRECT. Optimization
% Letters (2018). ISSN 1862-4472, 12 (7), 1699-1712,
% DOI: 10.1007/s11590-017-1228-4
%
% Constraint handling strategie taken from:
%--------------------------------------------------------------------------
% Stripinis, L., Paulavicius, R., Zilinskas, J.: "Penalty functions and
% two-step selection procedure based DIRECT-type algorithm for constrained
% global optimization".  Structural and Multidisciplinary Optimization,
% (2019). ISSN 1615-1488, DOI: 10.1007/s00158-018-2181-2
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[SS, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[MSS, CE, third, VAL] = Alocate(SS, VAL);

% Initialization step
[MV, MD, MSS, CE, VAL, minval, Xminval] = Initialization(Problem, MSS,...
    CE, VAL, third, SS);

while VAL.perror > SS.TOL                 % Main loop
%--------------------------------------------------------------------------
    % Selection of potential optimal hyper-rectangles step
    [Main, POH] = Find_poh(MV, MD);
    
    % Subdivide potential optimal hyper-rectangles
    [VAL, MSS, CE] = Calulcs(VAL, Problem, MSS, CE, third, POH,...
        Main, SS, 2);
    
    % Find minima values
    [MV, MD, minval, xatmin, VAL, Xminval] = Find_min(CE, MSS, VAL);
    
    % Update minima and check stopping conditions
    [VAL, SS] = Arewedone(minval, VAL, SS, MSS);
    
%--------------------------------------------------------------------------
end                                       % End of while

% Return value
minima      = minval;
if SS.G_nargout == 2
    xatmin    = abs(VAL.b - VAL.a).*Xminval(:, 1) + VAL.a;
elseif SS.G_nargout == 3
    xatmin    = abs(VAL.b - VAL.a).*Xminval(:, 1) + VAL.a;
    history = VAL.history(1:(VAL.itctr), 1:4);
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
    1e-8, 'globalmin', 0, 'globalxmin', 0);

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
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, third, VAL] = Alocate(SS, VAL)
%--------------------------------------------------------------------------
% Create Initial values
tic                                      % Mesure time
VAL.time    = 0;                         % initial time
VAL.fcount  = 1;                         % first fcnc counter
VAL.itctr   = 1;                         % initial iteration
VAL.perror  = 2*SS.TOL;                  % initial perror
VAL.epsil   = 1;                         % initial espil
VAL.tet     = 10^(-6);
VAL.ep      = 0.0001;
VAL.CARD    = 10*(VAL.n^3); % infeas HR
VAL.CARDi   = 1000*(VAL.n^3);
VAL.STAG    = 0;                         % iteration stagnates?
CE          = zeros(1, 10*VAL.n*SS.MAXdeep);% collum counter

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
    'L', zeros(VAL.n, 1), 'cc', zeros(1), 'ff', zeros(1));

third       = zeros(1, 10*VAL.n*SS.MAXdeep);      % delta values
third(1)    = 1/3;                       % first delta
for i = 2:10*VAL.n*SS.MAXdeep                     % all delta
    third(i)     = (1/3)*third(i - 1);
end

if SS.G_nargout == 3
    VAL.history  = zeros(SS.MAXits, 4);    % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [MV, MD, MSS, CE, VAL, minval, Xminval] = Initialization(...
    Problem, MSS, CE, VAL, third, SS)
%--------------------------------------------------------------------------
% Create Initial values
[MV, MD] = deal(ones(3, 1));                          % first fake value
MSS(1).L(:, 1) = zeros(VAL.n, 1);                      % Lengths
MSS(1).C(:, 1) = ones(VAL.n, 1)/2;                     % Center point
MSS(1).E(1) = 1;                                       % Index
point = abs(VAL.b - VAL.a).*MSS(1).C(:, 1) + VAL.a;    % Real point
[minval, MSS(1).F(1)] = deal(feval(Problem.f, point)); % f(x) eval.
MSS(1).cc(1) = CallC(Problem, point, SS);              % g(x) eval.
MSS(1).ff(1) = ~isequal(MSS(1).cc(1), 0);
CE(1) = 1; 

% Updtae iteration parameters
VAL.xatmin     = MSS(1).C(:, 1);
xatmin         = MSS(1).C(:, 1);
Xminval        = MSS(1).C(:, 1);

% Do we have feasible point?
if MSS.ff(1) ~= 0
    id = 0;
    fprintf('Phase II: searching feasible point:'); disp(' ');
    [MV, MD] = Phase(CE, MSS);
    while id == 0
%--------------------------------------------------------------------------
        % Selection of potential optimal hyper-rectangles step
        [Main, POH] = Find_poh(MV, MD);
        
        % Subdivide potential optimalhyper-rectangles
        [VAL, MSS, CE] = Calulcs(VAL, Problem, MSS, CE, third,...
            POH, Main, SS, 1);
        
        [MV, MD, minval] = Phase(CE, MSS);
        
        if minval  == 0
            id = 1;
            [MV, MD, minval, xatmin, VAL, Xminval] = Find_min(CE, MSS, VAL);
            fprintf('Violation of constraints: %15.10f    fn evals: %8i\n    f_min: %15.10f',...
                0, VAL.fcount, minval);
        else
            fprintf('Violation of constraints: %15.10f    fn evals: %8i\n',...
                minval, VAL.fcount);
        end
%--------------------------------------------------------------------------
    end
    fprintf('Phase I: Improve feasible solution:'); disp(' ');
end
VAL.xatmin   = xatmin;

% Check stop condition if global minima is known
VAL = Arewedone(minval, VAL, SS, MSS);

if SS.G_nargout == 3                     % Store History
    VAL.history(VAL.itctr, 1) = VAL.itctr;
    VAL.history(VAL.itctr, 2) = VAL.fcount;
    VAL.history(VAL.itctr, 3) = minval;
    VAL.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
return
%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [VAL, MSS, CE] = Calulcs(VAL, Problem, MSS, CE, third,...
    POH, Main, SS, P)
%--------------------------------------------------------------------------
[MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, Main, VAL, SS, P);

for i = size(POH, 2):-1:1                
    if ~isempty(POH{i})         
        if (CE(i) - size(POH{i}, 2)) == 0
            if find(CE ~= 0, 1, 'first') == i
                [MSS(i).E, MSS(i).L, MSS(i).F, MSS(i).C, MSS(i).cc,...
                    MSS(i).ff] = deal([]);
            end
        else
            C = setdiff(1:CE(i), POH{i});
            MSS(i).E(1:length(C)) = MSS(i).E(C);
            MSS(i).cc(1:length(C)) = MSS(i).cc(C);
            MSS(i).ff(1:length(C)) = MSS(i).ff(C);
            MSS(i).L(:, 1:length(C)) = MSS(i).L(:, C);
            MSS(i).F(1:length(C)) = MSS(i).F(C);
            MSS(i).C(:, 1:length(C)) = MSS(i).C(:, C);
        end
        CE(i) = CE(i) - size(POH{i}, 2);
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
%--------------------------------------------------------------------------
function [VAL, SS] = Arewedone(minval, VAL, SS, MSS)
%--------------------------------------------------------------------------
VAL.time = toc;

if SS.showITS == 1                % Show iteration stats
    fprintf(...
        'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, minval, VAL.time, VAL.fcount);
end

if SS.TESTflag == 1               % Check for stop condition
    if SS.globalMIN ~= 0            % Calculate error if globalmin known
        VAL.perror = 100*(minval - SS.globalMIN)/abs(SS.globalMIN);
    else
        VAL.perror = 100*minval;
    end
    if VAL.perror < SS.TOL
        fprintf('Minima was found with Tolerance: %4i', SS.TOL);
        VAL.perror = -10;
    end
else
    VAL.perror = 10;
end

if VAL.itctr >= SS.MAXits                % Have we exceeded the maxits?
    disp('Exceeded max iterations. Increase maxits');
    VAL.perror = -10;
end

if VAL.fcount > SS.MAXevals              % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

% Have we exceeded max deep?
if SS.MAXdeep <= max(MSS(end).L(:, 1)) + 1
    disp('Exceeded Max depth. Increse maxdeep'); VAL.perror = -1;
end

if SS.G_nargout == 3                     % Store History
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
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, Main, VAL,...
    SS, P)
%--------------------------------------------------------------------------
for i = 1:size(Main, 2)
    [DD, TT, mdx, VAL] = Evaluations(Problem, i, MSS, third, VAL,...
        Main, SS, P);
    for h = 1:TT
        II              = mdx + h;                % SET index
        if II > length(CE)
            CE(II) = 0;
        end
        if CE(II) == 0
            [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, CE(II - 1)*2));
            [MSS(II).F, MSS(II).E] = deal(zeros(1, CE(II - 1)*2));
            [MSS(II).cc, MSS(II).ff] = deal(zeros(1, CE(II - 1)*2));
        end
        IL              = CE(II) + 1;             % Left index
        IR              = CE(II) + 2;             % Right index
        CE(II)          = IR;                     % Colum index
        if CE(II) > size(MSS(II).F, 2)
            MSS(II).L     = [MSS(II).L, zeros(VAL.n, CE(II - 1)*2)];
            MSS(II).F     = [MSS(II).F, zeros(1, CE(II - 1)*2)];
            MSS(II).E     = [MSS(II).E, zeros(1, CE(II - 1)*2)];
            MSS(II).C     = [MSS(II).C, zeros(VAL.n, CE(II - 1)*2)];
            MSS(II).cc    = [MSS(II).cc, zeros(1, CE(II - 1)*2)];
            MSS(II).ff    = [MSS(II).ff, zeros(1, CE(II - 1)*2)];
        end
        MSS(II).F(IL)   = DD.L(DD.lsx(h));       % Left f(x)
        MSS(II).F(IR)   = DD.R(DD.lsx(h));       % Right f(x)
        MSS(II).E(IL)   = DD.eL(h);              % Left fcn counter
        MSS(II).E(IR)   = DD.eR(h);              % Right fcn counter
        MSS(II).L(:, IL) = DD.lL(:, h);          % Left lenght
        MSS(II).L(:, IR) = DD.lL(:, h);          % Right lenght
        MSS(II).C(:, IL) = DD.cL(:, DD.lsx(h));  % Left x
        MSS(II).C(:, IR) = DD.cR(:, DD.lsx(h));  % Right x
        MSS(II).cc(IL)  = DD.gL(DD.lsx(h));      % Left g(x)
        MSS(II).cc(IR)  = DD.gR(DD.lsx(h));      % Right g(x)
        MSS(II).ff(IL)  = DD.fL(DD.lsx(h));      % Left feasible flag
        MSS(II).ff(IR)  = DD.fR(DD.lsx(h));      % Right feasible flag
    end
    II                = mdx + TT;                  % SET index
    CE(II)            = CE(II) + 1;                % Colum index
    MSS(II).F(CE(II))     = DD.O;                  % Center f(x)
    MSS(II).E(CE(II))     = DD.eO;                 % Center fcn counter
    MSS(II).C(:, CE(II))  = DD.cO;                 % Center x
    MSS(II).L(:, CE(II))  = DD.LO;                 % Center lenght
    MSS(II).cc(CE(II))    = DD.gO;                 % Center g(x)
    MSS(II).ff(CE(II))    = DD.iO;                 % Center feasible flag
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function [Main, PH] = Find_poh(MM, DD)
%--------------------------------------------------------------------------
[s_a, s_b] = deal(0);
fc_min = MM(1, :);
dist_min = DD(1, :);
[POH_a, POH_b, PH] = deal(cell(1, size(fc_min, 2)));
[ss_a, ss_b] = deal(zeros(1, size(fc_min, 2)));
[index_a, index_b] = deal(size(fc_min, 2));
% Find index set of potential optimal hyper-rectangles
while index_a ~= 0
    [m_m, index_a] = min(fc_min(1:index_a));
    if ~isnan(m_m)
        s_a = s_a + 1;
        ss_a(s_a) = index_a;
        POH_a{index_a} = MM(2, index_a);
    end
    index_a = index_a - 1;
end
ss_a  = ss_a(1:s_a);

while index_b ~= 0
    [m_m, index_b] = min(dist_min(1:index_b));
    if ~isnan(m_m)
        s_b = s_b + 1;
        ss_b(s_b) = index_b;
        POH_b{index_b} = DD(2, index_b);
    end
    index_b = index_b - 1;
end
ss_b  = ss_b(1:s_b);

for i = 1:size(fc_min, 2)
   PH{i} = union(POH_a{i}, POH_b{i}); 
end

ss_b(ismember(DD(3, ss_b), intersect(DD(3, ss_b), MM(3, ss_a)))) = [];
Main = sortrows([[ss_a, ss_b]; MM(2, ss_a), DD(2, ss_b);...
    MM(3, ss_a), DD(3, ss_b)].', 3).';
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, S, mdx, VAL] = Evaluations(Const, i, MSS, third, VAL,...
    Main, SS, P)
%--------------------------------------------------------------------------
% Create and allocating empty sets
mdx   = Main(1, i);
A.cO  = MSS(mdx).C(:, Main(2, i));
A.LO  = MSS(mdx).L(:, Main(2, i));
A.eO  = MSS(mdx).E(Main(2, i));
A.O   = MSS(mdx).F(Main(2, i));
A.gO  = MSS(mdx).cc(Main(2, i));
A.iO  = MSS(mdx).ff(Main(2, i));
DIVIS = find(A.LO == min(A.LO));
DELTA = third(min(A.LO) + 1);
S     = length(DIVIS);
A.cL  = A.cO*ones(1, S);
A.cR  = A.cL;
A.lL  = A.LO*ones(1, S);
[A.L, A.R, A.gL, A.gR, A.eL, A.eR, A.fL, A.fR] =...
    deal(zeros(1, S));

for g = 1:S
    % left side
    A.cL(DIVIS(g), g) = A.cL(DIVIS(g), g) - DELTA;
    point = abs(VAL.b - VAL.a).*A.cL(:,g) + VAL.a;
    A.L(g) = feval(Const.f, point); A.gL(g) = CallC(Const, point, SS);
    A.fL(g) =  ~isequal(A.gL(g), 0);
    % right side
    A.cR(DIVIS(g), g) = A.cR(DIVIS(g), g) + DELTA;
    point = abs(VAL.b - VAL.a).*A.cR(:,g) + VAL.a;
    A.R(g) = feval(Const.f, point); A.gR(g) = CallC(Const, point, SS);
    A.fR(g) =  ~isequal(A.gR(g), 0);
end
if P == 2
    [~, A.lsx] = sort([min(A.L ,  A.R )' DIVIS], 1);
else
    [~, A.lsx] = sort([min(A.gL, A.gR)' DIVIS], 1);
end

for g = 1:S
    A.lL(DIVIS(A.lsx(1:g, 1)), g) = A.lL(DIVIS(A.lsx(1:g, 1)), g) + 1;
    A.eL(g) = VAL.fcount + 1;
    A.eR(g) = VAL.fcount + 2;
    VAL.fcount  = A.eR(g);
end
A.LO = A.lL(:, S);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  CallConstraints
% Purpose    :  Evaluate Constraints at pointed specified
%--------------------------------------------------------------------------
function ret_value = CallC(Problem, point, OPTI)
%--------------------------------------------------------------------------
ret_value = 0;
if isfield(Problem, 'constraint')
    [con_g, con_h] = feval(Problem.constraint, point);
    % inequality
    for i = 1:length(con_g)
        if con_g(i) > OPTI.ept
            ret_value = ret_value + con_g(i);
        end
    end
    % equality
    for i = 1:length(con_h)
        if abs(con_h(i)) > OPTI.ept
            ret_value = ret_value + abs(con_h(i));
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_min
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function [MV, MD, minval, xatmin, VAL, Xminval] = Find_min(CE, MSS, VAL)
%--------------------------------------------------------------------------
[MV, MD] = deal(nan(3, size(MSS, 2)));
ON  = 0;

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        TM = find(MSS(i).ff(1:CE(i)) == 0);
        if ~isempty(TM)
            TM = TM(MSS(i).E(TM) == max(MSS(i).E(TM(MSS(i).F(TM)...
            == min(MSS(i).F(TM))))));
            MV(1, i) = MSS(i).F(TM);
            MV(2, i) = TM;
            MV(3, i) = MSS(i).E(TM);
        end
    end
end
[minval, Least] = min(MV(1, :));
Xminval = MSS(Least).C(:,MV(2, Least));

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        TM = find(MSS(i).cc(1:CE(i)) > 0 & MSS(i).cc(1:CE(i)) < VAL.epsil);
        ON = ON + size(TM, 2);
        [MSS(i).ff(TM), MSS(i).cc(TM)] = deal(0);
        Temp = MSS(i).F(1:CE(i)) + MSS(i).cc(1:CE(i)) +...
            abs(MSS(i).F(1:CE(i)) - minval).*MSS(i).ff(1:CE(i));
        TMx = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(Temp == min(Temp))));
        MV(1, i) = Temp(TMx);
        MV(2, i) = TMx;
        MV(3, i) = MSS(i).E(TMx);
    end
end

if VAL.STAG == 10 && VAL.CARD <= VAL.CARDi
    VAL.STAG    = 0;
    VAL.epsil   = 1;
    VAL.CARD    = VAL.CARD * 10;
    VAL.tet     = VAL.tet/100;
elseif (ON == 0)       && (VAL.epsil * 3 <= 10)
    VAL.epsil   = VAL.epsil * 3;
elseif (ON > VAL.CARD) && (VAL.epsil / 3 >= VAL.ep)
    VAL.epsil   = VAL.epsil / 3;
elseif (ON > VAL.CARD) && (VAL.epsil / 3 <= VAL.ep)
    VAL.epsil = 0;
end

Least = find(MV(3, :) == max(MV(3, MV(1, :) == min(MV(1,:)))));
xatmin = MSS(Least).C(:, MV(2, Least));

if VAL.epsil == 0
    LEMDA = abs(((abs(VAL.b - VAL.a).*xatmin + VAL.a))...
            - (abs(VAL.b - VAL.a).*VAL.xatmin + VAL.a));
    if sum(LEMDA.^2)^0.5 < VAL.tet
        VAL.STAG  = VAL.STAG + 1;
    else
        VAL.STAG  = 0;
    end
end
VAL.xatmin = xatmin;

for i = 1:size(MSS, 2)                                     % DISTANCES
    if CE(i) ~= 0
        D_eul = sum((xatmin - MSS(i).C(:, 1:CE(i))).^2, 1).^0.5;
        TM = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(D_eul == min(D_eul))));
        MD(1, i) = D_eul(TM);
        MD(2, i) = TM;
        MD(3, i) = MSS(i).E(TM);
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Phase
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function [MV, MD, minval, xatmin] = Phase(CE, MSS)
%--------------------------------------------------------------------------
[MV, MD] = deal(nan(3, size(MSS, 2)));

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        TM = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(MSS(i).cc(1:CE(i))...
            == min(MSS(i).cc(1:CE(i))))));
        MV(1, i) = MSS(i).cc(TM);
        MV(2, i) = TM;
        MV(3, i) = MSS(i).E(TM);
    end
end

minval = min(MV(1, :));
Least = find(MV(1, :) == minval, 1, 'last');
xatmin = MSS(Least).C(:,MV(2, Least));

for i = 1:size(MSS, 2)                                     % DISTANCES
    if CE(i) ~= 0
        D_eul = sum((xatmin - MSS(i).C(:, 1:CE(i))).^2, 1).^0.5;
        TM = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(D_eul == min(D_eul))));
        MD(1, i) = D_eul(TM);
        MD(2, i) = TM;
        MD(3, i) = MSS(i).E(TM);
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