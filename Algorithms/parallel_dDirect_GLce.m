function [minima, xatmin, history] = parallel_dDirect_GLce(Problem,...
    opts, bounds)
%--------------------------------------------------------------------------
% Function   : parallel_dDirect_GLce
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 09/29/2019
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = parallel_dDirect_GLce(Problem, opts, bounds)
%       parallel  - Parallel implementation
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
%                 opts.cores     = number of computational threads
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
%
% Parallel scheme taken from:
%--------------------------------------------------------------------------
% Stripinis, L., Zilinskas, J., Casado, G. L., Paulavicius, R.:
% "Accelerating DIRECT-GLce algorithm for constrained global optimization
% through dynamic data structures and parallelization ".  Applied
% Mathematics and Computation, (2020).
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[SS, VAL, Problem] = Options(opts, nargout, Problem, bounds);

% Execute code in parallel on workers
spmd (SS.workers)
%--------------------------------------------------------------------------
    for i = 1:10000
        if labindex == 1
            labSend([], 2:SS.workers)
        else
            labReceive(1);
        end
    end
    % Create necessary data accessible to all workers
    [MSS, CE, third, VAL, Q] = Alocate(SS, VAL);
    
    if labindex == 1                    % master?
        [MSS, CE, Fmin, DATA, VAL, SS, XXMin] = Begin(Problem, MSS,...
            CE, VAL, SS);
    else                                % Receive information from master
        [DATA, ~] = labReceive(1);
    end
    
    while VAL.perror > SS.TOL           % Main loop
        [MSS, CE, Q] = Calulcs(VAL, Problem, DATA, MSS, CE, third, SS, Q);
        if labindex == 1                % Gather information from workers
            [DATA, Fmin, VAL, SS, XXMin] = MasterLAB(VAL, SS, MSS, Q);
        else                            % Receive information from master
            [DATA, ~,TAG] = labReceive(1);
            if TAG == 1, VAL.perror = -1; end
        end
    end                                  % End of while
%--------------------------------------------------------------------------
end                                     % End of SPMD block
minima  = Fmin{1};                      % Return value
xatmin  = XXMin{1};                     % Return point
TT      = SS{1};
VAL     = VAL{1};
history = TT.history(1:(VAL.itctr), 1:4); % Return history
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
    'testflag', 0, 'tol', 0.01, 'showits', 1, 'dimension', 1,...
    'globalmin', 0, 'globalxmin', 0, 'ept', 1e-8);

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

workers = gcp();
OPTI.workers = workers.NumWorkers;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Alocate
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, third, VAL, Q] = Alocate(SS, VAL)
%--------------------------------------------------------------------------
% Create Initial values
tic                                     % Mesure time
VAL.perror = 2*SS.TOL;                  % initial perror
VAL.aloc = 10*VAL.n*SS.MAXdeep;
CE = zeros(1, VAL.aloc);        % collum counter

% alociate MAIN sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
    'L', zeros(VAL.n, 1), 'cons', zeros(1), 'ff', zeros(1));
third = zeros(1, VAL.aloc);           % delta values
third(1) = 1/3;                         % first delta
for i = 2:VAL.aloc                    % all delta
    third(i) = (1/3)*third(i - 1);
end
Q = {ones(3, 1), [], [], [], ones(3, 1), [], [], ones(3, 1)};
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Begin
% Purpose   : Create necessary data for master
%--------------------------------------------------------------------------
function [MSS, CE, Fmin, DATA, VAL, SS, XXMin] = Begin(O, MSS, CE, VAL, SS)
%--------------------------------------------------------------------------
% Create and allocating empty sets
SS.history = zeros(SS.MAXits, 4);                      % allocating history

% Create Initial values
[VAL.fcn, VAL.itctr, MSS(1).E(1), CE(1)] = deal(1);
[VAL.time, VAL.STAG] = deal(0);                

% First evaluation
MSS(1).L(:, 1) = zeros(VAL.n, 1);                       % Lengths
[MSS(1).C(:, 1), VAL.xatmin, XXMin] = deal(ones(VAL.n, 1)/2);

Xmin = abs(VAL.b - VAL.a).*MSS(1).C(:, 1) + VAL.a;      % Real point
[Fmin, MSS(1).F(1)] = deal(feval(O.f, Xmin));           % f(x) eval.
MSS(1).cons(1) = CallConstraints(O, Xmin, SS.ept);      % g(x) eval.
MSS(1).ff(1) = ~isequal(MSS(1).cons(1), 0);

VAL.ep      = 10^(-4);
VAL.CARD    = 10*((size(VAL.a(:, 1), 1)))^3;            % POWER limit
VAL.PHASE   = isequal(MSS(1).cons(1), 0);               % PHASE I - II
VAL.CARDi   = 1000*((size(VAL.a(:, 1), 1)))^3; 
VAL.tet     = 10^(-6);
VAL.epsil   = 9;

% Divide values for workers
[MSV, MSD]  = deal(ones(4, 1));   % first fake value
[KK, TT, I, DEL, VAL.fcn] = Decision(1, 1, VAL, MSD, MSV, SS);

% work for master
DATA = {KK, TT, I, size(I, 1), DEL, VAL.xatmin, VAL.PHASE, Fmin, 9};

% Store History
if SS.G_nargout == 3
    SS.history(VAL.itctr, 1) = VAL.itctr;
    SS.history(VAL.itctr, 2) = VAL.fcn;
    SS.history(VAL.itctr, 3) = Fmin;
    SS.history(VAL.itctr, 4) = VAL.time;
end
%--------------------------------------------------------------------------
labSend(DATA, 2:SS.workers);                    % Send jobs to workers
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : MasterLAB
% Purpose   : last calculs of master
%--------------------------------------------------------------------------
function [DATA, fmin, VAL, SS, Xreturn] = MasterLAB(VAL, SS, MSS, Q)
%--------------------------------------------------------------------------
% collect information
[count, LB, CC, FC, EE, ES, MM, MN]  = deal(ones(1, SS.workers));
[VR, IR, ER, VF, IF, EF, DV, DI, DE] = deal(nan(SS.workers, VAL.aloc));
[XF, XR] = deal(zeros(VAL.n, SS.workers));
PW = 0;
for i = 1:SS.workers
    if i ~= labindex
        [Q, L] = labReceive('any');     % Receive info from workers
        LB(i) = L;                      % labindex
    end
    count(i)          = size(Q{1}, 2);
    PW                = PW + Q{7};            % power
    ES(i)             = Q{2}(3);
    FC(i)             = Q{6}(1);
    EE(i)             = Q{6}(2);
    CC(i)             = Q{4};
    MM(i)             = Q{2}(1);              % Values
    MN(i)             = Q{2}(2);
    VR(i, 1:count(i)) = Q{1}(1, :);           % Values
    IR(i, 1:count(i)) = Q{1}(2, :);           % indexes
    ER(i, 1:count(i)) = Q{1}(3, :);           % fcn
    VF(i, 1:count(i)) = Q{5}(1, :);           % Values
    IF(i, 1:count(i)) = Q{5}(2, :);           % indexes
    EF(i, 1:count(i)) = Q{5}(3, :);           % fcn
    XR(:, i)          = Q{3}(:, 1);           % x
    XF(:, i)          = Q{3}(:, 2);           % x
    DV(i, 1:count(i)) = Q{8}(1, :);
    DI(i, 1:count(i)) = Q{8}(2, :);
    DE(i, 1:count(i)) = Q{8}(3, :);
end
zz = max(count);

[MSV, MSC, MSD] = deal(nan(4, zz));
for i = 1:zz                        % Find min values
    mm1 = min(DV(:, i));
    mm2 = min(VR(:, i));
    if ~isnan(mm1)
        TM = find(DE(:, i) == max(DE(DV(:, i) == mm1, i)));
        MSD(1, i) = DV(TM, i);
        MSD(2, i) = DI(TM, i);
        MSD(3, i) = LB(TM);
        MSD(4, i) = DE(TM, i);
    end
    if ~isnan(mm2)    
        TM = find(ER(:, i) == max(ER(VR(:, i) == mm2, i)));
        MSV(1, i) = VR(TM, i);
        MSV(2, i) = IR(TM, i);
        MSV(3, i) = LB(TM);
        MSV(4, i) = ER(TM, i);
    end
end

% Create list KET of potentially optimal hyper-rectangles
POH_L = find_rectangles(MSD(1, :));
%--------------------------------------------------------------------------
if VAL.PHASE == 0
    MINI   = CC;
else
    MINI   = MM;
    if VAL.epsil ~= 0
        for i = 1:zz
            mm = min(VF(:, i));
            if ~isnan(mm)
                TM = find(EF(:, i) == max(EF(VF(:, i) == mm, i)));
                MSC(1, i) = VF(TM, i);
                MSC(2, i) = IF(TM, i);
                MSC(3, i) = LB(TM);
                MSC(4, i) = EF(TM, i);
            end
        end
    end
    
    if VAL.STAG == 10 && VAL.CARD <= VAL.CARDi
        VAL.STAG    = 0;
        VAL.epsil   = 1;
        VAL.CARD    = VAL.CARD * 10;
        VAL.tet     = VAL.tet/100;
    elseif (PW == 0)       && (VAL.epsil * 3 <= 10)
        VAL.epsil   = VAL.epsil * 3;
    elseif (PW > VAL.CARD) && (VAL.epsil / 3 >= VAL.ep)
        VAL.epsil   = VAL.epsil / 3;
    elseif (PW > VAL.CARD) && (VAL.epsil / 3 <= VAL.ep)
        VAL.epsil = 0;
    end
end

% Update f_min and x_min
fmin = min(MINI);                                 % f_min
Candidats = find(MINI == fmin);
[~, indexf] = max(EE(Candidats));
fminindex = (Candidats(indexf));

Xreturn = abs(VAL.b - VAL.a).*XR(:, fminindex) + VAL.a; % x_min
xmin = XR(:, fminindex);
constrain = CC(fminindex);                              % constrain

if VAL.PHASE == 0
    POH_G = find_rectangles(MSV(1, :));                 % POH
    [KK, TT, I, DEL, fcnc] = Decision(POH_G, POH_L, VAL, MSD, MSV, SS);
    DATA = {KK, TT, I, size(I, 1), DEL, XR(:, fminindex), VAL.PHASE, fmin, VAL.epsil};
    VAL.PHASE  = isequal(constrain, 0);                 % PHASE I - II
    if VAL.PHASE == 1
        fmin = FC(fminindex);                           % f_min
        DATA{8} = fmin;
        DATA{7} = 2;
    end
else
    if VAL.epsil ~= 0
        POH_G = find_rectangles(MSC(1, :));
        [KK, TT, I, DEL, fcnc] = Decision(POH_G, POH_L, VAL, MSD, MSC, SS);
        Candidats = find(MN == min(MN));
        [~, indexf] = max(ES(Candidats));
        fm = (Candidats(indexf));
        xmin = XF(:, fm);
    else
        POH_G = find_rectangles(MSV(1, :));
        [KK, TT, I, DEL, fcnc] = Decision(POH_G, POH_L, VAL, MSD, MSV, SS);
    end
        DATA = {KK, TT, I, size(I, 1), DEL, xmin, VAL.PHASE,...
        fmin, VAL.epsil};
end

if VAL.epsil == 0
    LEMDA = abs(((abs(VAL.b - VAL.a).*xmin + VAL.a))...
            - (abs(VAL.b - VAL.a).*VAL.xatmin + VAL.a));
    if sum(LEMDA.^2)^0.5 < VAL.tet
        VAL.STAG  = VAL.STAG + 1; else, VAL.STAG  = 0;
    end
    if VAL.STAG == 10, DATA{7} = 2; DATA{9} = 1; end
end
VAL.xatmin = xmin;

% Show iteration stats
if SS.showITS == 1
    VAL.time = toc;
    if constrain == 0
        fprintf('Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
            VAL.itctr,fmin,VAL.time,VAL.fcn);
    else
        fprintf('Iter: %4i   f_min: %15.10f*    time(s): %10.05f    fn evals: %8i\n',...
            VAL.itctr,fmin,VAL.time,VAL.fcn);
    end
end

% Check for stop condition
if SS.TESTflag == 1
    % Calculate error if globalmin known
    if VAL.PHASE == 1
        if SS.globalMIN ~= 0
            VAL.perror = 100*(fmin - SS.globalMIN)/abs(SS.globalMIN);
        else
            VAL.perror = 100*fmin;
        end
        if VAL.perror < SS.TOL
            fprintf('Minima was found with Tolerance: %4i', SS.TOL);
            VAL.perror = -1;
        end
    else
        VAL.perror = 10;
    end
else
    VAL.perror = 10;
end

% Have we exceeded the maxits?
if VAL.itctr >= SS.MAXits
    disp('Exceeded max iterations. Increase maxits');
    VAL.perror = -1;
end

% Have we exceeded the maxevals?
if VAL.fcn > SS.MAXevals
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -1;
end

% Have we exceeded max deep?
if  max(MSS(end).L(1))  > SS.MAXdeep
    disp('Exceeded Max depth. Increse maxdeep');
    VAL.perror = -1;
end

if VAL.perror == -1
    labSend(DATA, 2:SS.workers, 1);
else
    labSend(DATA, 2:SS.workers, 2);
end

% Store History
if SS.G_nargout == 3
    SS.history(VAL.itctr, 1) = VAL.itctr;
    SS.history(VAL.itctr, 2) = VAL.fcn;
    SS.history(VAL.itctr, 3) = fmin;
    SS.history(VAL.itctr, 4) = VAL.time;
end

% Update iteration number
VAL.fcn   = fcnc;
if VAL.perror > SS.TOL
    VAL.itctr = VAL.itctr + 1;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : find_rectangles
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function setas = find_rectangles(fc_min)
%--------------------------------------------------------------------------
s_i   = 0;
setas = zeros(1, size(fc_min, 2));
index = size(fc_min, 2);
% Find index set of potential optimal hyper-rectangles
while index ~= 0
    [m_m, index] = min(fc_min(1:index));
    if ~isnan(m_m)
        s_i = s_i + 1;
        setas(s_i) = index;
    end
    index = index - 1;
end
setas  = setas(1:s_i);
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [MSS, CE, Q] = Calulcs(VAL, Probl, DT, MSS, CE, third, SS, Q)
%--------------------------------------------------------------------------
% Create sets
if DT{9} ~= 0 && DT{7} == 1
    Q{1} = Q{5};
end
excess  = find(~cellfun(@isempty, DT{2})); % excess work on Workers
lack    = find(cellfun(@isempty, DT{2}));  % lack work on Workers
MAIN.FF = [];                                % create set MAIN

if ~isempty(excess)                  % Is any excess POH on Worker
    if ~isempty(DT{2}{labindex})
        % Create and allocating empty sets
        DO = size(DT{2}{labindex}, 2);
        MM = struct('F', zeros(1, DO), 'fef', zeros(1, DO), 'O',...
            zeros(1, DO), 'E', zeros(1, DO), 'C', zeros(VAL.n, DO), 'L',...
            zeros(VAL.n, DO), 'con', zeros(1, DO), 'ff', zeros(1, DO));
        for j = 1:DO                  % Write information which be send
            MM.O(j)    = DT{2}{labindex}(1, j);
            MM.con(j)  = MSS(MM.O(j)).cons(DT{2}{labindex}(2, j));
            MM.F(j)    = MSS(MM.O(j)).F(DT{2}{labindex}(2, j));
            MM.ff(j)   = MSS(MM.O(j)).ff(DT{2}{labindex}(2, j));
            MM.E(j)    = MSS(MM.O(j)).E(DT{2}{labindex}(2, j));
            MM.C(:, j) = MSS(MM.O(j)).C(:, DT{2}{labindex}(2, j));
            MM.L(:, j) = MSS(MM.O(j)).L(:, DT{2}{labindex}(2, j));
            MM.fef(j)  = DT{2}{labindex}(3, j);
        end
        labSend(MM, lack(1), 1);         % Send information to workers
    else
        % Create and allocating empty sets
        if labindex == lack(1)
            MM = struct('F', [], 'fef', [], 'O', [], 'E', [], 'C', [],...
                'L', [], 'con', [], 'ff', []);
            for j = 1:size(excess, 2)
                MOM    = labReceive(excess(j));
                MM.F   = [MM.F, MOM.F];
                MM.ff  = [MM.ff, MOM.ff];
                MM.O   = [MM.O, MOM.O];
                MM.C   = [MM.C, MOM.C];
                MM.L   = [MM.L, MOM.L];
                MM.E   = [MM.E, MOM.E];
                MM.fef = [MM.fef, MOM.fef];
                MM.con = [MM.con, MOM.con];
            end
        else
            MM = labReceive('any');
        end
        xc = DT{3}(labindex) - size(DT{1}{labindex}, 2);
        MAIN.FF  = MM.F(1:xc);      MM.F    = MM.F(xc + 1:end);
        MAIN.ff  = MM.ff(1:xc);     MM.ff   = MM.ff(xc + 1:end);
        MAIN.CC  = MM.C(:, 1:xc);   MM.C    = MM.C(:, xc + 1:end);
        MAIN.EE  = MM.E(1:xc);      MM.E    = MM.E(xc + 1:end);
        MAIN.LL  = MM.L(:, 1:xc);   MM.L    = MM.L(:, xc + 1:end);
        MAIN.OO  = MM.O(1:xc);      MM.O    = MM.O(xc + 1:end);
        MAIN.ISS = MM.fef(1:xc);    MM.fef  = MM.fef(xc + 1:end);
        MAIN.CON = MM.con(1:xc);    MM.con  = MM.con(xc + 1:end);
        if labindex ~= lack(end)
            labSend(MM, lack(find(lack == labindex) + 1));
        end
    end
end

% STORAGE of new rectangles
[MSS, CE] = STORAGE(Probl, MSS, CE, third, VAL, DT, MAIN, SS, Q{8}, Q{1});

% Find minima values
Q = find_min(CE, MSS, VAL.n, DT);
%--------------------------------------------------------------------------
labBarrier
if labindex ~= 1                      % Send info to master
    labSend(Q, 1);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS, CE] = STORAGE(O, MSS, CE, third, VAL, DATA, MAIN, SS, MD, MV)
%--------------------------------------------------------------------------
[DD, LIMIT] = Subdivision(O, DATA, MSS, third, VAL, MAIN, SS, DATA{7});

for i = 1:LIMIT                       % Calculations
    for h = 1:DD{i}.S
        II = DD{i}.mdx + h;                      % SET index
        if II > length(CE), CE(II) = 0; end
        if CE(II) == 0
            [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, 100));
            [MSS(II).F, MSS(II).E] = deal(zeros(1, 100));
            [MSS(II).cons, MSS(II).ff] = deal(nan(1, 100));
        end
        IL           = CE(II) + 1;        % Left index
        IR           = CE(II) + 2;        % Right index
        CE(II)       = IR;                % Colum index
        if CE(II) > size(MSS(II).F, 2)
            MSS(II).L    = [MSS(II).L, zeros(VAL.n, 100)];
            MSS(II).F    = [MSS(II).F, nan(1, 100)];
            MSS(II).cons = [MSS(II).cons, nan(1, 100)];
            MSS(II).ff   = [MSS(II).ff, zeros(1, 100)];
            MSS(II).E    = [MSS(II).E, zeros(1, 100)];
            MSS(II).C    = [MSS(II).C, nan(VAL.n, 100)];
        end
        MSS(II).F(IL) = DD{i}.L(DD{i}.lsx(h));        % Left f(x)
        MSS(II).F(IR) = DD{i}.R(DD{i}.lsx(h));        % Right f(x)
        MSS(II).ff(IL) = DD{i}.ffL(DD{i}.lsx(h));     % Left feas flag
        MSS(II).ff(IR) = DD{i}.ffR(DD{i}.lsx(h));     % Right feas flag
        MSS(II).cons(IL) = DD{i}.conl(DD{i}.lsx(h));  % Left sum g(x)
        MSS(II).cons(IR) = DD{i}.conr(DD{i}.lsx(h));  % Right sum g(x)
        MSS(II).E(IL) = DD{i}.eL(h);               % Left fcn counter
        MSS(II).E(IR) = DD{i}.eR(h);               % Right fcn counter
        MSS(II).L(:,IL) = DD{i}.lL(:, h);          % Left lenght
        MSS(II).L(:,IR) = DD{i}.lL(:, h);          % Right lenght
        MSS(II).C(:,IL) = DD{i}.cL(:, DD{i}.lsx(h));  % Left x
        MSS(II).C(:,IR) = DD{i}.cR(: ,DD{i}.lsx(h));  % Right x
    end
    II = DD{i}.mdx + DD{i}.S; 
    IO = CE(II) + 1;                 % Center index
    CE(II) = IO;                     % Colum index
    MSS(II).F(IO) = DD{i}.O;            % Center f(x)
    MSS(II).ff(IO) = DD{i}.ffo;         % Center feas flag
    MSS(II).E(IO) = DD{i}.eO;           % Center fcn counter
    MSS(II).cons(IO) = DD{i}.cono;       % Center sum g(x)
    MSS(II).C(:, IO) = DD{i}.cO;        % Center x
    MSS(II).L(:, IO) = DD{i}.lL(:, DD{i}.S); % Center lenght
end

D_O = DATA{5}{labindex};                % Indexes which can be replaced

for i = 1:size(D_O{1}, 2)                                                      
    I = D_O{1}(i);
    if size(find(D_O{2} == I), 2) ~= 0 && MD(2, I) == CE(I)                          
        MD(2,I) = MV(2,I);
    end

    if CE(I) == 0
        [MSS(I).E, MSS(I).L, MSS(I).C, MSS(I).F, MSS(I).ff, MSS(I).cons] = deal([]);
    else
        MSS(I).E([MV(2, I), CE(I)])     = MSS(I).E([CE(I), MV(2, I)]);                  
        MSS(I).L(:, [MV(2, I), CE(I)])  = MSS(I).L(:, [CE(I), MV(2, I)]);                
        MSS(I).F([MV(2, I), CE(I)])     = MSS(I).F([CE(I), MV(2, I)]);                   
        MSS(I).C(:, [MV(2, I), CE(I)])  = MSS(I).C(:, [CE(I), MV(2, I)]);   
        MSS(I).ff([MV(2, I), CE(I)])    = MSS(I).ff([CE(I), MV(2, I)]);
        MSS(I).cons([MV(2, I), CE(I)])  = MSS(I).cons([CE(I), MV(2, I)]);
    end
    CE(I) = CE(I) - 1;                                                    
end
for i = 1:size(D_O{2}, 2)                                                  
    I = D_O{2}(i);                                                                                                         
    if CE(I) == 0
        [MSS(I).E, MSS(I).L, MSS(I).C, MSS(I).F, MSS(I).ff, MSS(I).cons] = deal([]);
    else
        MSS(I).E([MD(2, I), CE(I)]) = MSS(I).E([CE(I), MD(2, I)]);                   
        MSS(I).L(:, [MD(2, I), CE(I)]) = MSS(I).L(:, [CE(I), MD(2, I)]);                  
        MSS(I).F([MD(2, I), CE(I)]) = MSS(I).F([CE(I), MD(2, I)]);                    
        MSS(I).C(:, [MD(2, I), CE(I)]) = MSS(I).C(:, [CE(I), MD(2, I)]);   
        MSS(I).ff([MD(2, I), CE(I)]) = MSS(I).ff([CE(I), MD(2, I)]);
        MSS(I).cons([MD(2, I), CE(I)]) = MSS(I).cons([CE(I), MD(2, I)]);
    end
    CE(I) = CE(I) - 1;                                                   
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Subdivision
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, LIMIT] = Subdivision(Problem, DATA, MSS, third, VAL, MAIN, SS, P)
%--------------------------------------------------------------------------
LIMIT_w = size(DATA{1}{labindex},2);
LIMIT_s = size(MAIN.FF, 2);
LIMIT = LIMIT_s + LIMIT_w;
A = cell(1, LIMIT);
fc = zeros(1, LIMIT);
LL = zeros(VAL.n, LIMIT);

for i = 1:LIMIT_w
    A{i}.mdx = DATA{1}{labindex}(1, i);
    A{i}.O   = MSS(A{i}.mdx).F(DATA{1}{labindex}(2, i));
    A{i}.eO  = MSS(A{i}.mdx).E(DATA{1}{labindex}(2, i));
    A{i}.cO  = MSS(A{i}.mdx).C(:, DATA{1}{labindex}(2, i));
    A{i}.cono = MSS(A{i}.mdx).cons(:, DATA{1}{labindex}(2, i));
    A{i}.ffo  = MSS(A{i}.mdx).ff(:, DATA{1}{labindex}(2, i));
    LL(:, i) = MSS(A{i}.mdx).L(:, DATA{1}{labindex}(2, i));
    fc(i)    = DATA{1}{labindex}(3, i);
end

for i = 1:LIMIT_s
    j = LIMIT_w + i;
    A{j}.mdx = MAIN.OO(i);
    A{j}.O = MAIN.FF(i);
    A{j}.eO = MAIN.EE(i);
    A{j}.cO = MAIN.CC(:, i);
    A{j}.cono = MAIN.CON(i);
    A{j}.ffo = MAIN.ff(i);
    LL(:, j) = MAIN.LL(:, i);
    fc(j) = MAIN.ISS(i);
end

for i = 1:LIMIT
    DIVIS = find(LL(:, i) == min(LL(:, i)));
    DELTA = third(min(LL(:, i)) + 1);
    A{i}.S = length(DIVIS);
    [A{i}.L, A{i}.R, A{i}.eL, A{i}.eR, A{i}.conl, A{i}.conr,...
        A{i}.ffL, A{i}.ffr] = deal(zeros(1, A{i}.S));
    A{i}.lL = LL(:, i)*ones(1, A{i}.S);
    A{i}.cL = A{i}.cO*ones(1, A{i}.S); A{i}.cR = A{i}.cL;
    
    for g = 1:A{i}.S
        A{i}.cL(DIVIS(g), g) = A{i}.cL(DIVIS(g), g) - DELTA;
        Point = abs(VAL.b - VAL.a).*A{i}.cL(:, g) + VAL.a;
        A{i}.L(g) = feval(Problem.f, Point);
        A{i}.conl(g) = CallConstraints(Problem, Point, SS.ept);
        if A{i}.conl(g) == 0, A{i}.ffL(g) = 0; else, A{i}.ffL(g) = 1; end
        
        A{i}.cR(DIVIS(g), g) = A{i}.cR(DIVIS(g), g) + DELTA;
        Point = abs(VAL.b - VAL.a).*A{i}.cR(:, g) + VAL.a;
        A{i}.R(g) = feval(Problem.f, Point);
        A{i}.conr(g) = CallConstraints(Problem, Point, SS.ept);
        if A{i}.conr(g) == 0, A{i}.ffR(g) = 0; else, A{i}.ffR(g) = 1; end
    end
    if P ~= 0
        [~, A{i}.lsx] = sort([min(A{i}.L, A{i}.R)' DIVIS], 1);
    else
        [~, A{i}.lsx] = sort([min(A{i}.conl, A{i}.conr)' DIVIS], 1);
    end
    
    for g = 1:A{i}.S
        A{i}.lL(DIVIS(A{i}.lsx(1:g, 1)), g) =...
                                   A{i}.lL(DIVIS(A{i}.lsx(1:g, 1)), g) + 1;
        A{i}.eL(g) = fc(i) + 1;
        A{i}.eR(g) = fc(i) + 2;
        fc(i) = fc(i) + 2;
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  find_min
% Purpose    :  Find min value
%--------------------------------------------------------------------------
function Q = find_min(CE, MSS, n, DT)
%--------------------------------------------------------------------------
% Create empty set
Q = {nan(3, size(MSS, 2)), nan(1, 3), zeros(n, 2), nan,...
    nan(3, size(MSS, 2)), nan(2, 1), 0, nan(3, size(MSS, 2))};
TT = nan(3, size(MSS, 2));

if size(MSS, 2) ~= 1
    for i = 1:size(MSS, 2)                                % DISTANCES
        if CE(i) ~= 0
            MSS(i).D = abs(MSS(i).F(1:CE(i)) - DT{8}).*MSS(i).ff(1:CE(i));
            D_eul = sum((DT{6}(:, 1) - MSS(i).C(:, 1:CE(i))).^2, 1).^0.5;
            TT(2, i) = find(MSS(i).E ==...
                           max(MSS(i).E(D_eul == min(D_eul))), 1, 'first');
            TT(1, i) = D_eul(TT(2, i));
            TT(3, i) = MSS(i).E(TT(2, i));
        end
    end
    POH = flip(find_rectangles(TT(1, :)));    % Find a local set of POH
    Q{8}(:, POH) = TT(:, POH);
    
    if DT{7} ~= 0                               % phase          
        if DT{9} ~= 0                           % epsil
            for i=1:size(MSS, 2)
                if CE(i) ~= 0
                    temp = MSS(i).D + MSS(i).cons(1:CE(i)) +...
                                                        MSS(i).F(1:CE(i));
                    TT(2, i) = find(MSS(i).E ==...
                             max(MSS(i).E(temp == min(temp))), 1, 'first');            
                    TT(1, i) = temp(TT(2, i));
                    TT(3, i) = MSS(i).E(TT(2, i));
                end
            end
            POH = flip(find_rectangles(TT(1, :))); % Find a local set of POH
            Q{1}(:, POH) = TT(:, POH);
            for i = 1:size(MSS, 2)
                if CE(i) ~= 0
                    temp = MSS(i).F(1:CE(i)) + (MSS(i).D +...
                            MSS(i).cons(1:CE(i))).*(MSS(i).cons(1:CE(i))...
                            >= DT{9});
                    TT(2, i) = find(MSS(i).E ==...
                             max(MSS(i).E(temp == min(temp))), 1, 'first');
                    TT(1, i) = temp(TT(2, i));
                    TT(3, i) = MSS(i).E(TT(2, i));
                    Q{7} = Q{7} + nnz((MSS(i).cons(1:CE(i)) <=...
                                        DT{9} & MSS(i).cons(1:CE(i)) > 0));
                end
            end
            POH = flip(find_rectangles(TT(1, :))); % Find a local set of POH
            Q{5}(:, POH) = TT(:, POH);
            
            [Q{2}(2), fmin_false] = min(Q{5}(1, :));
            Q{2}(3)    = MSS(fmin_false).E(:, Q{5}(2, fmin_false));
            Q{3}(:, 2) = MSS(fmin_false).C(:, Q{5}(2, fmin_false));
        else                                    % epsil
            for i = 1:size(MSS, 2)
                if CE(i) ~= 0
                    temp = MSS(i).D + MSS(i).F(1:CE(i)) +...
                          MSS(i).cons(1:CE(i)) + 0.0005*MSS(i).ff(1:CE(i));
                    TT(2, i) = find(MSS(i).E ==...
                             max(MSS(i).E(temp == min(temp))), 1, 'first');            
                    TT(1, i) = temp(TT(2, i));
                    TT(3, i) = MSS(i).E(TT(2, i));
                end
            end
            POH = flip(find_rectangles(TT(1, :))); % Find a local set of POH
            Q{1}(:, POH) = TT(:, POH);
        end
        
    else                                        % phase
        for i = 1:size(MSS, 2)
            if CE(i) ~= 0      
                TT(2, i) = find(MSS(i).E ==...
                    max(MSS(i).E(MSS(i).cons(1:CE(i)) ==...
                    min(MSS(i).cons(1:CE(i))))), 1, 'first');            
                TT(1, i) = MSS(i).cons(TT(2, i));
                TT(3, i) = MSS(i).E(TT(2, i));
            end
        end
        POH = flip(find_rectangles(TT(1, :)));    % Find a local set of POH
        Q{1}(:, POH) = TT(:, POH);

    end
    [Q{2}(1), index] = min(Q{1}(1, :));
    Q{6}(1) = MSS(index).F(:, Q{1}(2, index));
    Q{6}(2) = MSS(index).E(:, Q{1}(2, index));
    Q{3}(:, 1) = MSS(index).C(:, Q{1}(2, index));
    Q{4} = MSS(index).cons(:, Q{1}(2, index));
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  CallConstraints
% Purpose    :  Evaluate Constraints at pointed specified
%--------------------------------------------------------------------------
function ret_value = CallConstraints(Problem, point, ept)
%--------------------------------------------------------------------------
ret_value = 0;
if isfield(Problem, 'constraint')
    [con_g, con_h] = feval(Problem.constraint, point);
    % inequality
    for i = 1:length(con_g)
        if con_g(i) > ept
            ret_value = ret_value + con_g(i);
        end
    end
    % equality
    for i = 1:length(con_h)
        if abs(con_h(i)) > ept
            ret_value = ret_value + abs(con_h(i));
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Decision
% Purpose   : Decide which workers whould delete theirs rectnagles
%--------------------------------------------------------------------------
function [KK, TT, I, HH, fnc] = Decision(K, L, VAL, MSD, MSV, SS)
%--------------------------------------------------------------------------
L(ismember(MSD(4, L), intersect(MSV(4, K), MSD(4, L)))) = [];
Main = [K, L; MSV(2, K), MSD(2, L);...
            VAL.fcn, zeros(1, size([K, L], 2) - 1); MSV(4, K), MSD(4, L)];                
                         
I = flipud(split_scalar((size(Main, 2)), SS.workers));
[HH, KK, II, TT] = deal(cell(1, SS.workers));
I(size(I, 1) + 1:SS.workers) = 0;

for i = 2:size(Main, 2)
    [Main(3, i), VAL.fcn] = deal(VAL.fcn + (VAL.n -...
                                    mod(Main(1, i - 1) - 1, VAL.n))*2);
end
fnc = (VAL.fcn + (VAL.n - mod(Main(1, end) - 1, VAL.n))*2);

% Calculations
for i=1:SS.workers
    II{i} = find([MSV(3, K), MSD(3, L)] == i);
    HH{i}{1} = K(MSV(3, K) == i);
    HH{i}{2} = L(MSD(3, L) == i);
    if size(II{i}, 2) ~= 0
        if size(II{i}, 2) < (I(i)+1)
            KK{i} = Main(:, II{i});
        else
            TT{i} = Main(:, II{i}(I(i) + 1:size(II{i}, 2)));
            KK{i} = Main(:, II{i}(1:I(i)));
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Split scalar
% Purpose    :  Divides rectangle i that is passed in
%--------------------------------------------------------------------------
function [integerPerTask, numTasks] = split_scalar(intVal, numTasks)
%--------------------------------------------------------------------------
narginchk(2, 2);

if intVal > 0 && numTasks == 0
    error('pctexample:splitscalar:SplitScalarInvalidNumTasks', ...
        ['Number of tasks must be greater than 0 if the scalar is '...
        'greater than 0']);
end
if intVal < numTasks
    numTasks = intVal;
end

split                         = fix(intVal/numTasks);
remainder                     = intVal - numTasks*split;
integerPerTask                = zeros(numTasks, 1);
integerPerTask(:)             = split;
integerPerTask(1:remainder)   = integerPerTask(1:remainder) + 1;
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