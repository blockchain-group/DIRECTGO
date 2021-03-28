function [minima, xatmin, history] = dDirect_GLh(Problem, bounds, opts)
%--------------------------------------------------------------------------
% Function   : dDirect_GLh
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 06/15/2020
% Purpose    : DIRECT optimization algorithm for problems with various
%              constraints constraints
%--------------------------------------------------------------------------
% [minima, xatmin, history] = dDirect_GLh(Problem, bounds, opts)
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       G         - Enhancing the global search
%       L         - Enhancing the local search
%       h         - Hidden constraint handling
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f           = Objective function handle
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
%--------------------------------------------------------------------------

% Get options
SS = Options(opts, nargout);

% Alocate sets and create initial variables
[MSS, CE, third, VAL] = Alocate(bounds, SS);

% Initialization step
[MV, MD, MSS, CE, VAL, minval, xatmin] = Initialization(Problem, MSS,...
    CE, VAL, SS, third);

while VAL.perror > SS.TOL                 % Main loop
%--------------------------------------------------------------------------
    % Selection of potential optimal hyper-rectangles step
    [Main, POH] = Find_poh(MV, MD);

    % Subdivide potential optimalhyper-rectangles
    [VAL, MSS, CE, MV, minval, xatmin, MD] = Calulcs(VAL, Problem,...
        MSS, CE, third, POH, Main, minval, xatmin, 2, SS);
    
    % Update minima and check stopping conditions
    VAL = Arewedone(MSS, minval, VAL, SS, CE);
    
%--------------------------------------------------------------------------
end                                       % End of while

% Return value
minima      = minval;
if SS.G_nargout == 3
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
function SS = Options(opts, narg)
%--------------------------------------------------------------------------
% Determine option values
if nargin < 3 && (isempty(opts))
    opts = [];
end
getopts(opts, 'maxits', 1000, 'maxevals', 100000, 'maxdeep', 2000,...
    'testflag', 0, 'globalmin', 0, 'tol', 0.01, 'showits', 1, 'ept', 0);

SS.G_nargout = narg;     % output arguments
SS.MAXits    = maxits;   % maximum of iterations
SS.MAXevals  = maxevals; % maximum # of function evaluations
SS.MAXdeep   = maxdeep;  % maximum number of side divisions
SS.showITS   = showits;  % print iteration stat
SS.TOL       = tol;      % allowable relative error if f_reach is set
SS.TESTflag  = testflag; % terminate if within a relative tolerence of f_opt
SS.globalMIN = globalmin;% minimum value of function
SS.ept       = ept;      % tollerance for constraints
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, third, VAL] = Alocate(bounds, SS)
%--------------------------------------------------------------------------
% Create Initial values
tic                                         % Mesure time
VAL.a          = bounds(:, 1);              % left bound
VAL.b          = bounds(:, 2);              % right bound
VAL.n          = size(bounds, 1);           % dimension
VAL.time       = 0;                         % initial time
VAL.fcount     = 1;                         % first fcnc counter
VAL.itctr      = 1;                         % initial iteration
VAL.perror     = 2*SS.TOL;                  % initial perror
VAL.max        = sum((zeros(VAL.n, 1) - ones(VAL.n, 1)).^2, 1).^0.5;
CE             = zeros(1, 10*VAL.n*SS.MAXdeep);% collum counter

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
    'L', zeros(VAL.n, 1), 'X', ones(1));

third           = zeros(1, 10*VAL.n*SS.MAXdeep);     % delta values
third(1)        = 1/3;                               % first delta
for i = 2:10*VAL.n*SS.MAXdeep                        % all delta
    third(i)    = (1/3)*third(i - 1);
end
if SS.G_nargout == 3
    VAL.history = zeros(SS.MAXits, 4);        % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Initialization of the DIRECT
%--------------------------------------------------------------------------
function [MV, MD, MSS, CE, VAL, Fmin, Xmin] = Initialization(...
    Problem, MSS, CE, VAL, SS, third)
%--------------------------------------------------------------------------
% Create Initial values
[MV, MD] = deal(ones(3,1));                            % first fake value
MSS(1).L(:, 1) = zeros(VAL.n, 1);                      % Lengths
MSS(1).C(:, 1) = ones(VAL.n, 1)/2;                     % Center point
MSS(1).E(1) = 1;                                       % Index
MSS(1).diam = 1/2*norm((1/3*(ones(VAL.n, 1))).^MSS(1).L(:, 1));
Point = abs(VAL.b - VAL.a).*MSS(1).C(:, 1) + VAL.a;
Xmin = MSS(1).C(:, 1);
CE(1)         = 1;
MSS.X(1)      = CallC(Problem, Point, SS);
if MSS.X(1) == 0
    MSS.F(1)    = feval(Problem.f, Point);
else
    MSS.F(1)    = 10^9;
end
Fmin          = MSS(1).F(1);                           % initial minima

if MSS.X(1) == 1                       % initial midpoint feasible?
%--------------------------------------------------------------------------
    POH{1}      = 1;
    Main        = ones(3, 1);
    LH_index    = 1;
    fprintf('Phase II: searching feasible point');
    
    while LH_index ~= 0
        [VAL, MSS, CE, ~, Fmin, Xmin] = Calulcs(VAL, Problem,...
            MSS, CE, third, POH, Main, Fmin, Xmin, 1, SS);
        
        % find largest hyper-rectangles, with largest index
        POH = cell(1, size(MSS, 2));
        Main = zeros(3, 1);
        idx = find(CE ~= 0, 1, 'first');
        [~, II] = max((MSS(idx).E(1:CE(idx))));
        Main(1, 1) = idx;
        Main(2, 1) = II;
        Main(3, 1) = MSS(idx).E(II);
        POH{idx} = II;
        % do we have a feasible midpoint?
        if ~isempty(find([MSS.X] == 0, 1))
            LH_index = 0;
        end
        fprintf('fn evals: %8i\n', VAL.fcount);
    end
    [MV, MD, Fmin, Xmin, MSS] = Find_min(CE, MSS, VAL, 2);
%--------------------------------------------------------------------------
end

% Check stop condition if global minima is known
VAL = Arewedone(MSS, Fmin, VAL, SS, CE);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Arewedone
% Purpose   : Update minima value and check stopoing conditions
%--------------------------------------------------------------------------
function VAL = Arewedone(MSS, minval, VAL, SS, CE)
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
if SS.MAXdeep < max(MSS(find(CE ~= 0, 1, 'last')).L(:, 1))
    disp('Exceeded Max depth. Increse maxdeep'); VAL.perror = -1;
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
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [VAL, MSS, CE, MV, minval, xatmin, MD] = Calulcs(VAL,...
    Problem, MSS, CE, third, POH, Main, minval, xatmin, ph, OPTI)
%--------------------------------------------------------------------------
[MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, VAL, Main,...
    minval, xatmin, ph, OPTI);

for i = size(POH, 2):-1:1                
    if ~isempty(POH{i})         
        if (CE(i) - size(POH{i}, 2)) == 0
            if find(CE ~= 0, 1, 'first') == i
                [MSS(i).E, MSS(i).L, MSS(i).F, MSS(i).C] = deal([]);
            end
        else
            C = setdiff(1:CE(i), POH{i});
            MSS(i).E(1:length(C)) = MSS(i).E(C);
            MSS(i).L(:, 1:length(C)) = MSS(i).L(:, C);
            MSS(i).F(1:length(C)) = MSS(i).F(C);
            MSS(i).C(:, 1:length(C)) = MSS(i).C(:, C);
            MSS(i).X(:, 1:length(C)) = MSS(i).X(:, C);
        end
        CE(i) = CE(i) - size(POH{i}, 2);
    end
end

% Find minima values
[MV, MD, minval, xatmin, MSS] = Find_min(CE, MSS, VAL, ph, xatmin);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = STORAGE(Problem, MSS, CE, third, VAL,...
    Main, minval, xatmin, ph, OPTI)
%--------------------------------------------------------------------------
for i = 1:size(Main, 2)
    [DD, TT, mdx, VAL] = Evaluations(Problem, i, MSS, third, VAL,...
        Main, minval, xatmin, ph, OPTI);
    for h = 1:TT
        II              = mdx + h;               % SET index
        if CE(II) == 0
            MSS(II).diam = 1/2*norm((1/3*(ones(VAL.n, 1))).^DD.lL(:, h));
            [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, CE(II - 1)*2));
            [MSS(II).F, MSS(II).E, MSS(II).X] = deal(ones(1, CE(II - 1)*2));
        end
        IL              = CE(II) + 1;            % Left index
        IR              = CE(II) + 2;            % Right index
        CE(II)          = IR;                    % Colum index
        if CE(II) > size(MSS(II).F, 2)
            MSS(II).L      = [MSS(II).L, zeros(VAL.n, CE(II - 1)*2)];
            MSS(II).F      = [MSS(II).F, nan(1, CE(II - 1)*2)];
            MSS(II).E      = [MSS(II).E, zeros(1, CE(II - 1)*2)];
            MSS(II).C      = [MSS(II).C, zeros(VAL.n, CE(II - 1)*2)];
            MSS(II).X      = [MSS(II).X, nan(1, CE(II - 1)*2)];
        end
        MSS(II).F(IL)    = DD.L(DD.lsx(h, 1));   % Left f(x)
        MSS(II).F(IR)    = DD.R(DD.lsx(h, 1));   % Right f(x)
        MSS(II).X(IL)    = DD.Lx(DD.lsx(h, 1));  % Left flag
        MSS(II).X(IR)    = DD.Rx(DD.lsx(h, 1));  % Right flag
        MSS(II).E(IL)    = DD.eL(h);             % Left fcn counter
        MSS(II).E(IR)    = DD.eR(h);             % Right fcn counter
        MSS(II).L(:, IL) = DD.lL(:, h);          % Left lenght
        MSS(II).L(:, IR) = DD.lL(:, h);          % Right lenght
        MSS(II).C(:, IL) = DD.cL(:,DD.lsx(h, 1));% Left x
        MSS(II).C(:, IR) = DD.cR(:,DD.lsx(h, 1));% Right x
    end
    II                    = mdx + TT;          % SET index
    CE(II)                = CE(II) + 1;        % Colum index
    MSS(II).F(CE(II))     = DD.O;              % Center f(x)
    MSS(II).X(CE(II))     = DD.Ox;             % Center flag
    MSS(II).E(CE(II))     = DD.eO;             % Center fcn counter
    MSS(II).C(:, CE(II))  = DD.cO;             % Center x
    MSS(II).L(:, CE(II))  = DD.LO;             % Center lenght
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Evaluations
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, S, mdx, VAL] = Evaluations(Const, i, MSS, third, VAL,...
    Main, Fmin, Xmin, ph, OPTI)
%--------------------------------------------------------------------------
% Create and allocating empty sets
mdx = Main(1, i);
A.cO = MSS(mdx).C(:,Main(2, i));
A.LO = MSS(mdx).L(:,Main(2, i));
A.eO = MSS(mdx).E(Main(2, i));
A.O = MSS(mdx).F(Main(2, i));
A.Ox = MSS(mdx).X(Main(2, i));
DIVIS = find(A.LO == min(A.LO));
DELTA = third(min(A.LO) + 1);
S = length(DIVIS);
A.lL = A.LO*ones(1, S);
A.cL = A.cO*ones(1, S);
A.cR = A.cL;
[A.L, A.R, A.Lx, A.Rx, A.eL, A.eR] = deal(zeros(1, S));

for g = 1:S
    A.cL(DIVIS(g), g) = A.cL(DIVIS(g), g) - DELTA;
    om_point = abs(VAL.b - VAL.a).*A.cL(:, g) + VAL.a;
    A.Lx(g) = CallC(Const, om_point, OPTI);
    if A.Lx(g) == 0
        A.L(g) = feval(Const.f, om_point);
    else
        A.L(g) = 10^9;
    end
    if A.Lx(g) == 1 && ph == 2
        A.L(g) = Fmin + sum((Xmin - om_point).^2, 1).^0.5;
    end
    
    A.cR(DIVIS(g), g) = A.cR(DIVIS(g), g) + DELTA;
    om_point = abs(VAL.b - VAL.a).*A.cR(:, g) + VAL.a;
    A.Rx(g) = CallC(Const, om_point, OPTI);
    if A.Rx(g) == 0
        A.R(g)    = feval(Const.f, om_point);
    else
        A.R(g)    = 10^9;
    end
    if A.Rx(g) == 1 && ph == 2
        A.R(g) = Fmin + sum((Xmin - om_point).^2, 1).^0.5;
    end
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
% Function   :  CallConstraints
% Purpose    :  Evaluate Constraints at pointed specified
%--------------------------------------------------------------------------
function ff_value = CallC(Problem, point, OPTI)
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
if ret_value == 0
    ff_value = 0;
else
    ff_value = 1;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Find_min
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function [MV, MD, minval, xatmin, MSS] = Find_min(CE, MSS, VAL, ph, xatmin)
%--------------------------------------------------------------------------
[MV, MD] = deal(nan(3, size(MSS, 2)));

if ph == 1
    minval = 10^(9);
    maxval = minval + 1;
else
    MA = nan(4, size(MSS, 2));
    for i = 1:size(MSS, 2)
        if CE(i) ~= 0
            TX = find(MSS(i).X(1:CE(i)) == 0);
            if ~isempty(TX)
                TM = TX(MSS(i).F(TX) == min(MSS(i).F(TX)));
                if size(TM, 2) ~= 1
                    [~, II]  = max((MSS(i).E(TM)));
                    TM = TM(II);
                end
                MA(1, i) = MSS(i).F(TM);
                MA(2, i) = TM;
                MA(3, i) = MSS(i).E(TM);
                MA(4, i) = max(MSS(i).F(1:CE(i)));
            end
        end
    end
    
    [minval, index] = min(MA(1, :));
    maxval = max(MA(4, :));
    xatmin = MSS(index).C(:, MA(2, index));
end

for i = 1:size(MSS, 2)                                     % DISTANCES
    if CE(i) ~= 0
        D_eul = sum((xatmin - MSS(i).C(:, 1:CE(i))).^2, 1).^0.5;
        TM = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(D_eul == min(D_eul))));
        MD(1, i) = D_eul(TM);
        MD(2, i) = TM;
        MD(3, i) = MSS(i).E(TM);

        TX = find(MSS(i).X(1:CE(i)) == 1);
        MSS(i).F(TX) = minval + D_eul(TX);
        TT = abs(MSS(i).F(1:CE(i)) - minval  + 10^(-16))/abs(maxval...
            - minval  + 10^(-16));
        TT(TX) = D_eul(TX)/VAL.max;
        TT = TT*MSS(i).diam;

        TM = find(MSS(i).E(1:CE(i)) == max(MSS(i).E(TT == min(TT))));
        MV(1, i) = TT(TM);
        MV(2, i) = TM;
        MV(3, i) = MSS(i).E(TM);
    end
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