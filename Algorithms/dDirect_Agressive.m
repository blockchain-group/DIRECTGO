function [minima, xatmin, history] = dDirect_Agressive(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : dDirect_Agressive
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dDirect_Agressive(Problem, opts, bounds)
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       Agressive - aggressive selection of potential optimal rectangles
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
% Selection of potential optimal hyper-rectangles taken from:
%--------------------------------------------------------------------------
% Baker, C.A., Watson, L.T., Grossman, B., Mason, W.H., Haftka, R.T.
% "Parallel global aircraft con?guration design space exploration".
% In: A. Tentner (ed.) High Performance Computing Symposium 2000,
% pp. 54–66. Soc. for Computer Simulation Internat (2000)
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[SS, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[MSS, CE, third, VAL] = Alocate(SS, VAL);

% Initialization step
[MSS, CE, VAL, minval, xatmin, POH] = Initialization(Problem,...
    MSS, CE, VAL, SS);

while VAL.perror > SS.TOL                 % Main loop
    % Subdivide potential optimalhyper-rectangles
    [VAL, MSS, CE, minval, xatmin, POH] = Calulcs(VAL, Problem,...
        MSS, CE, third, POH);
    
    % Update minima and check stopping conditions
    [VAL, SS] = Arewedone(MSS, minval, VAL, SS);
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
    'testflag', 0, 'tol', 0.01, 'showits', 1, 'dimension', 1,...
    'globalmin', 0, 'globalxmin', 0);

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
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, third, VAL] = Alocate(SS, VAL)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
VAL.time       = 0;                         % initial time
VAL.fcount     = 1;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
CE             = zeros(1, VAL.n*SS.MAXdeep);% collum counter

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
    'L', zeros(VAL.n, 1));

third           = zeros(1, SS.MAXdeep);     % delta values
third(1)        = 1/3;                      % first delta
for i = 2:SS.MAXdeep                        % all delta
    third(i)    = (1/3)*third(i - 1);
end
if SS.G_nargout == 3
    VAL.history = zeros(SS.MAXits, 4);       % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [MSS, CE, VAL, Fmin, Xmin, POH] = Initialization(Problem,...
    MSS, CE, VAL, SS)
%--------------------------------------------------------------------------
% Create Initial values
MSS(1).L(:, 1) = zeros(VAL.n, 1);                 	% Lengths
MSS(1).Diam = 1/2*norm((1/3*(ones(VAL.n, 1))).^MSS(1).L(:, 1));
MSS(1).C(:, 1) = ones(VAL.n, 1)/2;               	% Center point
MSS(1).E(1) = 1;                                    % Index
MSS(1).F(1) = feval(Problem.f, abs(VAL.b - VAL.a).*(MSS(1).C(:, 1))...
    + VAL.a);
CE(1)         = 1;

[Fmin, Xmin, POH] = Find_min(CE, MSS);
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
function [VAL, MSS, CE, minval, xatmin, POH] = Calulcs(VAL,...
    Problem, MSS, CE, third, POH)
%--------------------------------------------------------------------------
[MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, VAL, POH);

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
[minval, xatmin, POH] = Find_min(CE, MSS);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, VAL, PH)
%--------------------------------------------------------------------------
for i = 1:size(MSS, 2)
    if ~isempty(PH{i})
        for j = 1:size(PH{i}, 2)
            [DD, TT, VAL] = CALCULS_WORKER(Problem, i, MSS, third,...
                VAL, PH{i}(j));
            for h = 1:TT
                II = i + h;      
                if II > size(CE, 2)
                  CE(II) = 0;    
                end
                if CE(II) == 0
                    MSS(II).Diam = 1/2*...
                        norm((1/3*(ones(VAL.n, 1))).^DD.lL(:, h));
                    [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, 100));
                    [MSS(II).F, MSS(II).E] = deal(zeros(1, 100));
                end
                IL = CE(II) + 1;             
                IR = CE(II) + 2;             
                CE(II) = IR;                    
                if CE(II) > size(MSS(II).F, 2)
                    MSS(II).L = [MSS(II).L, zeros(VAL.n, 100)];
                    MSS(II).F = [MSS(II).F, zeros(1, 100)];
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
function [A, S, VAL] = CALCULS_WORKER(Const, i, MSS, third, VAL, PH)
%--------------------------------------------------------------------------
% Create and allocating empty sets
A.cO = MSS(i).C(:, PH);
A.LO = MSS(i).L(:, PH);
A.eO = MSS(i).E(PH);
A.O = MSS(i).F(PH);
DIVIS = find(A.LO == min(A.LO));
DELTA = third(min(A.LO(:, 1)) + 1);
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
    A.lL(DIVIS(A.lsx(1:g, 1)), g) = A.lL(DIVIS(A.lsx(1:g, 1)), g) + 1;
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
function [VAL, SS] = Arewedone(MSS, minval, VAL, SS)
%--------------------------------------------------------------------------
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

% Have we exceeded max deep?
if SS.MAXdeep <= max(MSS(end).L(:, 1)) + 1
    disp('Exceeded Max depth. Increse maxdeep'); VAL.perror = -1;
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
% Function   :  Find_min
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function [minval, xatmin, POH] = Find_min(CE, MSS)
%--------------------------------------------------------------------------
MV = nan(1, size(MSS, 2));
POH = cell(1, size(MSS, 2));

for i = 1:size(MSS, 2)
    if CE(i) ~= 0
        TM = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(MSS(i).F(1:CE(i))...
            == min(MSS(i).F(1:CE(i))))));
        MV(1, i) = MSS(i).F(TM);
        POH{i} = TM;
    end
end

[minval, Least] = min(MV(1, :));
xatmin = MSS(Least).C(:, POH{Least});
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