function [ret_minval, final_xatmin, history] = parallel_dDirect_sym2...
    (Problem, bounds, opts)
%--------------------------------------------------------------------------
% Function   : parallel_dDirect_sym2
% Author 1   : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2   : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 04/19/2020
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = parallel_dDirect_sym2(Problem, bounds, opts)
%       parallel  - parallel version of the algorithm
%       d         - dynamic memory management in data structure
%       DIRECT    - DIRECT(DIvide a hyper-RECTangle) algorithm
%       sym       - DIRECT extension for symmetric functions
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
%                 opts.showits   = 1 print iteration status
%                 opts.ep        = global/local weight parameter
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
% Grbic, R., Nyarko, E.K., Scitovski, R. "A modifcation of the DIRECT
% method for Lipschitz global optimization for a symmetric function".
% J. Global Optim. 1�20 (2012). DOI 10.1007/s10898-012-0020-3
%--------------------------------------------------------------------------

% Get options
if nargin < 3
    opts = [];
end
SS = Options(opts, nargout);

spmd                              % Execute code in parallel on labs
%--------------------------------------------------------------------------
    % Alocate sets and create initial variables
    [MSS, CE, VAL] = Alocate(bounds, SS);
    
    if labindex == 1 
        [MSS, CE, DT, VAL, SS, xmin, fmin] = Initialization(Problem,...
            MSS, CE, VAL, SS);
    else        
        DT = labReceive(1);
    end
    while VAL.perror > SS.TOL       % Main loop
        [MSS, CE, DT] = subdivision(VAL, Problem, DT, MSS, CE, SS);
        if labindex == 1            % Gather information from workers
            [DT, fmin, VAL, SS, xmin] = Master(VAL, SS, DT, MSS);
        else                        % Receive information from master
            [DT, ~, TAG] = labReceive(1); 
            if TAG == 1
                VAL.perror = -1; 
            end
        end
    end                              % End of while
%--------------------------------------------------------------------------
end                                 % End of SPMD block
ret_minval      = fmin{1};          % Return value
final_xatmin    = xmin{1};       	% Return point
TT              = SS{1};
VAL             = VAL{1};
history         = TT.history(1:(VAL.itctr - 1), 1:4); % Return history
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
getopts(opts, 'maxits', 1000, 'maxevals', 100000, 'maxdeep', 1000,...
    'testflag', 0, 'globalmin', 0, 'tol', 0.01, 'showits', 1, 'ep', 1e-4);

SS.G_nargout = narg;     % output arguments
SS.MAXits    = maxits;   % maximum of iterations
SS.MAXevals  = maxevals; % maximum # of function evaluations
SS.MAXdeep   = maxdeep;  % maximum number of side divisions
SS.showITS   = showits;  % print iteration stat
SS.TOL       = tol;      % allowable relative error if f_reach is set
SS.TESTflag  = testflag; % terminate if within a relative tolerence of f_opt
SS.globalMIN = globalmin;% minimum value of function
SS.ep        = ep;       % global/local weight parameter

workers = gcp();
SS.workers = workers.NumWorkers;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : BEGIN
% Purpose   : Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [MSS, CE, VAL] = Alocate(bounds, SS)
%--------------------------------------------------------------------------
% Create Initial values
tic                                     % Mesure time
VAL.a = bounds(:, 1);                   % left bound
VAL.b = bounds(:, 2);                   % right bound
VAL.n = size(bounds, 1);                % dimension
VAL.count = 1;
VAL.perror = 2*SS.TOL;                  % initial perror
CE = zeros(1, VAL.n*SS.MAXdeep);        % collum counter

% Create and allocating empty sets
MSS = struct('F', zeros(1), 'E', zeros(1), 'C', zeros(VAL.n, 1),...
    'L', zeros(VAL.n, 1));
L = ones(VAL.n, 1)/2;
MSS(1).Diam = sqrt(sum(L.^2));
for i = 1:SS.MAXdeep
    k = VAL.n*(i - 1);
    for j = 1:VAL.n
        L(j) = L(j)/3;
        MSS(k + j + 1).Diam = round(sqrt(sum(L.^2)), 15);
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Initialization
% Purpose   : Create necessary data for master
%--------------------------------------------------------------------------
function [MSS, CE, DATA, VAL, SS, xmin, fmin] = Initialization(Problem,...
    MSS, CE, VAL, SS)
%--------------------------------------------------------------------------
% Create and allocating empty sets
SS.history = zeros(SS.MAXdeep, 4);              % allocating history
VAL.fcount = 1;                                 % first fcnc counter
VAL.itctr = 1;                                  % initial iteration
VAL.time = 0;                                   % initial time
CE(1) = 1;                                      % \# in column.

% First evaluation
MSS(1).L(:, 1) = ones(VAL.n, 1)/2;              % Lengths
MSS(1).C(:, 1) = ones(VAL.n,1)/2;               % Center point
MSS(1).E(1) = 1;                                % Index
[point, xmin] = deal(abs(VAL.b-VAL.a).*MSS(1).C(:, 1) + VAL.a);
[MSS(1).F(1), fmin] = deal(feval(Problem.f, point));                      

% Divide values for workers
[KK, TT, I, DEL, VAL.fcount] = DIVas(VAL, SS, {1; 1; 1; 1; 1});
% work for master
DATA = {KK, TT, I, size(I,1), DEL, MSS(1).F(1)};

%--------------------------------------------------------------------------
labSend(DATA, 2:SS.workers);               % Send jobs to workers
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : subdivision
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [MSS, CE, DATA] = subdivision(VAL, Problem, DT, MSS, CE, SS)
%--------------------------------------------------------------------------
% Create sets
excess = find(~cellfun(@isempty, DT{2}));       % excess work on Workers
lack = find(cellfun(@isempty, DT{2}));          % lack work on Workers
MAIN.FF = [];                                   % create set MAIN

if ~isempty(excess)                            % Is any excess POH on Worker
    if ~isempty(DT{2}{labindex})
        % Create and allocating empty sets
        DO = size(DT{2}{labindex}, 2);
        [MM.F, MM.O, MM.E, MM.fef] = deal(zeros(1, DO));
        [MM.C, MM.L] = deal(zeros(VAL.n, DO));
        for j = 1:DO
            % Write information which be send
            MM.F(j) = MSS(DT{2}{labindex}(1, j)).F(DT{2}{labindex}(2, j));
            MM.E(j) = MSS(DT{2}{labindex}(1, j)).E(DT{2}{labindex}(2, j));
            MM.C(:, j) = MSS(DT{2}{labindex}(1, j)).C(:, DT{2}{labindex}(2, j));
            MM.L(:, j) = MSS(DT{2}{labindex}(1, j)).L(:, DT{2}{labindex}(2, j));
            MM.fef(j) = DT{2}{labindex}(3, j);
            MM.O(j) = DT{2}{labindex}(1, j);
        end
        
        labSend(MM, lack(1));       	 % Send information to other labs
    else
        % Receive information
        if labindex == lack(1)
            [MO.F, MO.O, MO.E, MO.fef, MO.L, MO.C] = deal([]);
            for j = 1:size(excess, 2)

                MOM = labReceive('any');
                MO.F = [MO.F, MOM.F];
                MO.O = [MO.O, MOM.O];
                MO.C = [MO.C, MOM.C];
                MO.L = [MO.L, MOM.L];
                MO.E = [MO.E, MOM.E];
                MO.fef = [MO.fef, MOM.fef];
            end
        else
            MO = labReceive('any');
        end
        jj = find(lack == labindex);
        xc = DT{3}(lack(jj)) - size(DT{1}{lack(jj)}, 2);
        MAIN.FF = MO.F(1:xc); MO.F = MO.F(xc + 1:end);
        MAIN.CC = MO.C(:, 1:xc); MO.C = MO.C(:, xc + 1:end);
        MAIN.EE = MO.E(1:xc);   MO.E = MO.E(xc + 1:end);
        MAIN.LL = MO.L(:, 1:xc); MO.L = MO.L(:, xc + 1:end);
        MAIN.OO = MO.O(1:xc);   MO.O = MO.O(xc + 1:end);
        MAIN.ISS = MO.fef(1:xc); MO.fef = MO.fef(xc + 1:end);
        if size(lack, 2) ~= jj
            labSend(MO, lack(jj + 1));
        end
    end
end

% STORAGE of new rectangles
[MSS, CE] = STORAGE(Problem, MSS, CE, VAL, DT, MAIN);

% Find minima values
[xmin, fmin, hul] = find_min(CE, MSS, VAL.n, SS, DT{6});
DATA = {xmin, fmin, hul};
%--------------------------------------------------------------------------
labBarrier
if labindex ~= 1                        % Send info to master
    labSend(DATA, 1);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : STORAGE
% Purpose   : Store information wrom workers
%--------------------------------------------------------------------------
function [MSS,CE] = STORAGE(Problem, MSS, CE, VAL, DATA, MAIN)
%--------------------------------------------------------------------------
[DD, LIMIT] = Calulcs(Problem, DATA, MSS, VAL, MAIN);

for i = 1:LIMIT
    for h = 1:DD{i}.S
        II = DD{i}.mdx + h;                          
        if CE(II) == 0
            [MSS(II).L, MSS(II).C] = deal(zeros(VAL.n, 100));
            [MSS(II).F, MSS(II).E] = deal(zeros(1, 100));
        end
        if CE(II) > size(MSS(II).F,2)
            MSS(II).L = [MSS(II).L, zeros(VAL.n, 100)];
            MSS(II).F = [MSS(II).F, nan(1, 100)];
            MSS(II).E = [MSS(II).E, zeros(1, 100)];
            MSS(II).C = [MSS(II).C, zeros(VAL.n, 100)];
        end
        discard = symdirect(VAL, DD{i}.cL(:, DD{i}.lsx(h)), DD{i}.lL(:, h));
        if discard == false
            IL           = CE(II) + 1;
            CE(II)       = IL; 
            MSS(II).F(IL) = DD{i}.L(DD{i}.lsx(h));
            MSS(II).E(IL) = DD{i}.eL(h); 
            MSS(II).L(:, IL) = DD{i}.lL(:, h); 
            MSS(II).C(:, IL) = DD{i}.cL(:, DD{i}.lsx(h)); 
        end
        discard = symdirect(VAL, DD{i}.cR(:, DD{i}.lsx(h)), DD{i}.lL(:, h));
        if discard == false
            IR           = CE(II) + 1;
            CE(II)       = IR; 
            MSS(II).F(IR) = DD{i}.R(DD{i}.lsx(h));
            MSS(II).E(IR) = DD{i}.eR(h); 
            MSS(II).L(:, IR) = DD{i}.lL(:, h); 
            MSS(II).C(:, IR) = DD{i}.cR(:, DD{i}.lsx(h));
        end

    end
    discard = symdirect(VAL, DD{i}.cO, DD{i}.lL(:, DD{i}.S));
    if discard == false
        II = DD{i}.mdx + DD{i}.S;                       % SET index
        CE(II) = CE(II) + 1;                            % Colum index
        MSS(II).F(CE(II)) = DD{i}.O;                    % Center f(x)
        MSS(II).E(CE(II)) = DD{i}.eO;                   % Center fcn counter
        MSS(II).C(:, CE(II)) = DD{i}.cO;                % Center x
        MSS(II).L(:, CE(II)) = DD{i}.lL(:, DD{i}.S);    % Center lenght
    end
end

D_O = DATA{5}{1}{labindex};                % Indexes which can be replaced
D_Os = DATA{5}{2}{labindex};
remm = unique(D_O);

if isempty(remm) == 0
    for i = 1:size(remm, 2)
        I = remm(i);
        remov = D_Os(D_O == I);
        if (CE(I) - size(remov, 2)) == 0
            if find(CE ~= 0, 1, 'first') == I
                [MSS(I).E, MSS(I).L, MSS(I).F, MSS(I).C] = deal([]);
            end
        else
            C = setdiff(1:CE(I), remov);
            MSS(I).E(1:length(C)) = MSS(I).E(C);
            MSS(I).L(:, 1:length(C)) = MSS(I).L(:, C);
            MSS(I).F(1:length(C)) = MSS(I).F(C);
            MSS(I).C(:, 1:length(C)) = MSS(I).C(:, C);
        end
        CE(I) = CE(I) - size(remov, 2);
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Calulcs
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [A, LIMIT] = Calulcs(PR, DATA, MSS, VAL, MAIN)
%--------------------------------------------------------------------------
LIMIT_w = size(DATA{1}{labindex},2);
LIMIT_s = size(MAIN.FF, 2);
LIMIT = LIMIT_s + LIMIT_w;
A = cell(1, LIMIT);
fc = zeros(1, LIMIT);
LL = zeros(VAL.n, LIMIT);

for i = 1:LIMIT_w
    A{i}.mdx = DATA{1}{labindex}(1, i);
    A{i}.O = MSS(A{i}.mdx).F(DATA{1}{labindex}(2, i));
    A{i}.eO = MSS(A{i}.mdx).E(DATA{1}{labindex}(2, i));
    A{i}.cO = MSS(A{i}.mdx).C(:, DATA{1}{labindex}(2, i));
    LL(:, i) = MSS(A{i}.mdx).L(:, DATA{1}{labindex}(2, i));
    fc(i) = DATA{1}{labindex}(3, i);
end

for i = 1:LIMIT_s
    j = LIMIT_w + i;
    A{j}.mdx = MAIN.OO(i);
    A{j}.O = MAIN.FF(i);
    A{j}.eO = MAIN.EE(i);
    A{j}.cO = MAIN.CC(:, i);
    LL(:, j) = MAIN.LL(:, i);
    fc(j) = MAIN.ISS(i);
end

for i = 1:LIMIT
    max_L  = max(LL(:, i));
    DIVIS = find(LL(:, i) == max_L);
    DELTA = 2*max_L/3;
    A{i}.S = length(DIVIS);
    [A{i}.L, A{i}.R, A{i}.eL, A{i}.eR] = deal(zeros(1, A{i}.S));
    A{i}.lL = LL(:, i)*ones(1, A{i}.S);
    A{i}.cL = A{i}.cO*ones(1, A{i}.S); A{i}.cR = A{i}.cL;
    
    for g = 1:A{i}.S
        A{i}.cL(DIVIS(g), g) = A{i}.cL(DIVIS(g), g) - DELTA;
        A{i}.L(g) = feval(PR.f, abs(VAL.b - VAL.a).*A{i}.cL(:, g) + VAL.a);
        A{i}.cR(DIVIS(g), g) = A{i}.cR(DIVIS(g), g) + DELTA;
        A{i}.R(g) = feval(PR.f, abs(VAL.b - VAL.a).*A{i}.cR(:, g) + VAL.a);
    end
    [~, A{i}.lsx] = sort([min(A{i}.L, A{i}.R)' DIVIS], 1);
    for g = 1:A{i}.S
        A{i}.lL(DIVIS(A{i}.lsx(1:g, 1)), g) = DELTA/2;
                                  
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
function [xmin, fmin, hul] = find_min(CE, MSS, n, SS, minval)
%--------------------------------------------------------------------------
if isempty(find(CE ~= 0, 1)) == 0
    % Create sets   
    length_of_MSS = find(CE ~= 0, 1, 'last');
    MV = nan(2, length_of_MSS);
    nonsize  = [MSS.Diam];
    hul = cell(5, length_of_MSS);
    hs = nan(2, length_of_MSS);

    for i = 1:length_of_MSS
        if CE(i) ~= 0
            [MV(1, i), MV(2, i)] = min(MSS(i).F(1:CE(i)));
            ties = find(abs(MSS(i).F(1:CE(i)) -  MV(1, i)) <= 1E-12);
            if ~isempty(ties)
                hul{1, i} = ties;
                hul{2, i} = nonsize(i)*ones(1, size(ties, 2));
                hul{3, i} = MSS(i).F(ties);
                hul{4, i} = i*ones(1, size(ties, 2));
                hul{5, i} = labindex*ones(1, size(ties, 2));
            end
        end
    end

    [fmin, index] = min(MV(1, :));
    xmin = MSS(index).C(:,MV(2, index));
    minval = min([fmin, minval]);
    for i = 1:length_of_MSS
        if CE(i) ~= 0
            [hs(1, i), hs(2, i)] = min((MSS(i).F(1:CE(i)) - minval +...
                max(SS.ep*abs(minval), 1E-8))/nonsize(i));
    
        end
    end
else
    xmin = zeros(n, 1);
    fmin = nan(1, 2);
    hul = cell(5, 1);
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function   :  Split scalar
% Purpose    :  Divides rectangle i that is passed in
%--------------------------------------------------------------------------
function integerPerTask = split_scalar(intVal, numTasks)
%--------------------------------------------------------------------------
if intVal < numTasks
    numTasks = intVal;
end

split                         = fix(intVal / numTasks);
remainder                     = intVal - numTasks * split;
integerPerTask                = zeros(numTasks, 1);
integerPerTask(:)             = split;
integerPerTask(1:remainder)   = integerPerTask(1:remainder) + 1;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Master
% Purpose   : last calculs of master
%--------------------------------------------------------------------------
function [DATA, fmin, VAL, SS, xmin] = Master(VAL, SS, D, MSS)
%--------------------------------------------------------------------------
% Block any action until all workers reach this point
[count, LB, MN]  = deal(ones(1, SS.workers));
XR = zeros(VAL.n,SS.workers);
[Md, Mz, Mx, My, Mt] = deal(cell(1, VAL.count + VAL.n));
for i=1:SS.workers
    if i ~= labindex
        [D, L] = labReceive('any');     % Receive info from workers
        LB(i) = L;                      % labindex
    end
    count(i) = size(D{3}, 2);
    XR(:, i) = D{1}(:, 1);              % x
    MN(i) = D{2}(1);                    % minval
    for j = 1:count(i)
        Md{1, j} = [cell2mat(Md(1, j)), cell2mat(D{3}(1, j))];
        Mz{1, j} = [cell2mat(Mz(1, j)), cell2mat(D{3}(2, j))];
        Mx{1, j} = [cell2mat(Mx(1, j)), cell2mat(D{3}(3, j))];
        My{1, j} = [cell2mat(My(1, j)), cell2mat(D{3}(4, j))];
        Mt{1, j} = [cell2mat(Mt(1, j)), cell2mat(D{3}(5, j))];
    end
end

% Update f_min and x_min
[fmin, fminindex] = min(MN);                         % f_min
xmin = abs(VAL.b-VAL.a).*XR(:, fminindex) + VAL.a;
VAL.count = max(count);
MVV = nan(1, VAL.count);

for i = 1:VAL.count                     % Find min values
    if ~cellfun(@isempty, Mx(i))
        MVV(i) = min(cell2mat(Mx(i)));
        ties = find(abs(cell2mat(Mx(i)) -  MVV(i)) <= 1E-12);
        Md{i} = Md{i}(ties);
        Mz{i} = Mz{i}(ties);
        Mx{i} = Mx{i}(ties);
        My{i} = My{i}(ties);
        Mt{i} = Mt{i}(ties);
    end
end

nonsize  = [MSS.Diam];
hul = cell(5, VAL.count);
hs = nan(2, VAL.count);

for i = 1:VAL.count
    if ~isempty(Md{i})
        [hs(1, i), hs(2, i)] = min((Mx{i} - fmin +...
            max(SS.ep*abs(fmin), 1E-8))/nonsize(i));
        
    end
end
i_min = find(hs(1, :) == min(hs(1, :)));

POH = Find_po([Md; Mz; Mx; My], i_min, nonsize, VAL.count, MVV,...
    Mx{i_min}(hs(2, i_min)));
      
for i = 1:length(POH)
    hul{1, POH(i)} = Md{1, POH(i)};
    hul{2, POH(i)} = Mz{1, POH(i)};
    hul{3, POH(i)} = Mx{1, POH(i)};
    hul{4, POH(i)} = My{1, POH(i)};
    hul{5, POH(i)} = Mt{1, POH(i)};
end

[KK, TT, I, DEL, fcnc] = DIVas(VAL, SS, hul);
DATA = {KK, TT, I, size(I, 1), DEL, fmin};

% Show iteration stats
if SS.showITS == 1
    VAL.time = toc;
    fprintf('Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr, fmin, VAL.time, VAL.fcount);
end
% Check for stop condition
if SS.TESTflag == 1
    % Calculate error if globalmin known
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
% Have we exceeded the maxits?
if VAL.itctr >= SS.MAXits
    disp('Exceeded max iterations. Increase maxits'); VAL.perror = -1;
end
% Have we exceeded the maxevals?
if VAL.fcount > SS.MAXevals
    disp('Exceeded max fcn evals. Increase maxevals'); VAL.perror = -1;
end
% Have we exceeded max deep?
if SS.MAXdeep <= VAL.count + VAL.n
    disp('Exceeded Max depth. Increse maxdeep'); VAL.perror = -1;
end
if VAL.perror == -1
    labSend(DATA, 2:SS.workers,1);
else
    labSend(DATA, 2:SS.workers,2);
end
% Store History
if SS.G_nargout == 3
    SS.history(VAL.itctr,1) = VAL.itctr;
    SS.history(VAL.itctr,2) = VAL.fcount;
    SS.history(VAL.itctr,3) = fmin;
    SS.history(VAL.itctr,4) = VAL.time;
end

% Update iteration number
VAL.fcount = fcnc;
VAL.itctr = VAL.itctr + 1;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Decision
% Purpose   : Decide which workers whould delete theirs rectnagles
%--------------------------------------------------------------------------
function [KK, TT, I, HH, fnc] = DIVas(VAL, SS, HUL)
%--------------------------------------------------------------------------
K = cell2mat(HUL(4, :));
index_i = cell2mat(HUL(1, :));
MSV = cell2mat(HUL(5, :));
Main  = [K; cell2mat(HUL(1, :)); VAL.fcount, zeros(1, size(K, 2) - 1)];

I = flipud(split_scalar((size(Main, 2)), SS.workers));
[HH, KK, II, TT] = deal(cell(1, SS.workers));
I(size(I, 1) + 1:SS.workers) = 0;

for i = 2:size(Main, 2)
    [Main(3, i), VAL.fcount] = deal(VAL.fcount + (VAL.n -...
                                    mod(Main(1, i - 1) - 1, VAL.n))*2);
end
fnc = (VAL.fcount + (VAL.n - mod(Main(1, end) - 1, VAL.n))*2);

% Calculations
for i=1:SS.workers
    ii = MSV == i;
    II{i} = find(ii);
    HH{1}{i} = K(ii);
    HH{2}{i} = index_i(ii);
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
% Function  : Find_poh
% Purpose   : Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function S = Find_po(hulls, i_min, d, index, d_min, b1)
% Find all rects on hub
jj    = find(~cellfun(@isempty, hulls(1, :)), 1, 'first');
S_1   = cell(1, index);
S_1(1:i_min)   = hulls(1, 1:i_min);
Su    = cell(4, index);

if i_min - jj > 1
    a1 = d(i_min);
    a2 = d(jj);
    b2 = d_min(jj);
    % The line is defined by: y = slope*x + const
    slope = (b2 - b1)/(a2 - a1);
    const = b1 - slope*a1;
    for i = 1:i_min
        if ~isempty(hulls{3, i})
            TT = find(hulls{3, i} <= slope*d(i) + const + 1E-12);
            if ~isempty(TT)
                Su{1, i} = hulls{1, i}(TT);
                Su{2, i} = hulls{2, i}(TT);
                Su{3, i} = hulls{3, i}(TT);
                Su{4, i} = hulls{4, i}(TT);
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
    S = find(~cellfun(@isempty, S));
else
    S = find(~cellfun(@isempty, S_1));
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
m = length(x);
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
START = 1;
v = START;
w = length(x);
flag = 0;
% Index vector for points in convex hull
h = (1:length(x))';
while (next(v ,m)~=START) || (flag==0)
    if next(v, m) == w
        flag = 1;
    end
    a = v;
    b = next(v, m);
    c = next(next(v, m) ,m);
    if det([x(a) y(a) 1 ; x(b) y(b) 1 ; x(c) y(c) 1]) >= 0
        leftturn = 1;
    else
        leftturn = 0;
    end
    if leftturn
        v = next(v, m);
    else
        j = next(v, m);
        x = [x(1:j - 1); x(j + 1:m)];
        y = [y(1:j - 1); y(j + 1:m)];
        h = [h(1:j - 1); h(j + 1:m)];
        m = m - 1;
        w = w - 1;
        v = pred(v, m);
    end
end
xy = [x, y];
hh = size(xy, 1);
Z = cell(1, hh);
for i = 1:hh
    Z{i} = find(xy(i, 1) == xyAR(:, 1) & xy(i, 2) == xyAR(:, 2))';
end
h = cell2mat(Z);
return

function i = next(v, m)
if v == m
    i = 1;
else
    i = v + 1;
end
return

function i = pred(v, m)
if v == 1
    i = m;
else
    i = v - 1;
end
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
    equal = 0;
    if check == VAL.n - 1
        for j = 1:VAL.n - 1
            if V(j, i) - V(j + 1, i) < 1e-12
                equal = equal + 1;
            end
        end
        if equal ~= VAL.n - 1
            discard = false; % Cannot discard rectangle
            break;
        end
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