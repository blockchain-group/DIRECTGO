function [minima, xatmin, history] = HALRECT(Problem, opts, bounds)
%--------------------------------------------------------------------------
% Function   : HALRECT
% Written by : Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Written by : Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
% Created on : 04/30/2022
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% [minima, xatmin, history] = HALRECT(Problem, opts, bounds)
%       HALRECT - Hyper-rectangular partitioning based on 1-Dimensional 
%                Bisection and objective function evaluations at the center
%       G      - Enhancing the global search
%       L      - Enhancing the local search
%       IO     - Improved original selection scheme
%       IA     - Improved aggressive selection scheme
%
% Input parameters:
%       Problem - Structure containing problem
%                 Problem.f       = Objective function handle
%
%       opts    - MATLAB structure which contains options.
%                 opts.maxevals   = max. number of function evals
%                 opts.maxits     = max. number of iterations
%                 opts.maxdeep    = max. number of rect. divisions
%                 opts.globalmin  = globalmin (if known)
%                 opts.globalxmin = globalxmin (if known)
%                 opts.testflag   = 1 if globalmin known, 0 otherwise
%                 opts.dimension  = problem dimension
%                 opts.showits    = 1 print iteration status
%                 opts.ep         = global/local weight parameter
%                 opts.tol        = tolerance for termination if
%                                  testflag = 1
%                 opts.poh        = Choose the the selection of POH scheme:
%                           'GL' - Two-step-based Pareto selection (default) 
%                           'IA' - Improved Aggressive selection strategy 
%                           'IO' - Improved Original selection strategy
%                 opts.equation   = Selection is performed using one of 
%                                   the four posible models: '          ', 
%                           '22b', '22c', '22d'(default)
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
% Stripinis, L., Paulavicius, R., Zilinskas, J.: Improved scheme for
% selection of potentially optimal hyperrectangles in DIRECT. Optimization
% Letters (2018). ISSN 1862-4472, 12 (7), 1699-1712,
% DOI: 10.1007/s11590-017-1228-4
%--------------------------------------------------------------------------
% Baker, C.A., Watson, L.T., Grossman, B., Mason, W.H., Haftka, R.T.
% "Parallel global aircraft con?guration design space exploration".
% In: A. Tentner (ed.) High Performance Computing Symposium 2000,
% pp. 54–66. Soc. for Computer Simulation Internat (2000)
%--------------------------------------------------------------------------
% Jones, D.R.: Direct global optimization algorithm. In: Floudas, C.A., 
% Pardalos, P.M. (Eds.) Encyclopedia of Optimization, pp. 431–440. 
% Springer, Boston (2001). https://doi.org/10.1007/0-306-48332-7_93
%--------------------------------------------------------------------------
if nargin == 2, bounds = []; end
if nargin == 1, bounds = []; opts = []; end

% Get options
[OPTI, VAL] = Options(opts, nargout, Problem, bounds);

% Alocate sets and create initial variables
[VAL, MSS] = Alocate(OPTI, VAL);

% Initialization step
[OPTI, VAL, MSS] = Inits(VAL, OPTI, Problem, MSS);

while VAL.perror > OPTI.TOL                                 % Main loop
    % Selection of potential optimal hyper-rectangles step
    switch OPTI.poh
        case 'IO'
            aa = find(MSS.DD(1:VAL.I) ~= -1);
            POH = aa(Find_po(MSS.FE(aa), VAL.Fmin, OPTI.ep, MSS.DD(aa)));
        case 'IA'
            aa = find(MSS.DD(1:VAL.I) ~= -1);
            POH = aa(Aggressive(MSS.FE(aa), MSS.DD(aa), 1:length(aa), VAL, MSS));
        otherwise
            POH = Pareto(VAL, MSS);
    end

    [MSS, VAL] = Subdivisionas(VAL, Problem, MSS, POH);

    % Update minima and check stopping conditions
    [VAL, MSS] = Arewedone(OPTI, VAL, MSS);
end                                                         % End of while

% Return value
minima      = VAL.Fmin;
if OPTI.G_nargout == 2
    xatmin    = (abs(VAL.b - VAL.a)).*VAL.Xmin(:, 1) + VAL.a;
elseif OPTI.G_nargout == 3
    xatmin    = (abs(VAL.b - VAL.a)).*VAL.Xmin(:, 1) + VAL.a;
    history   = VAL.history(1:VAL.itctr, 1:4);
end

%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% AUXILIARY FUNCTION BLOCK
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Function: Aggressive
% Purpose:  Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function boxes = Aggressive(fc, szes, indexs, VAL, MSS)
%--------------------------------------------------------------------------
C = unique(szes);
if length(C) > VAL.n*100
    C = C((end - (VAL.n*100 - 1)):end);
end

SMS   = zeros(3, size(C, 2));
m_set = [fc; szes; indexs];
for i = 1: size(C, 2)                             % reduce m_set
    MB          = m_set(:, m_set(2, :) == C(i));
    [MV, MI]    = min(MB(1, :));
    TM          = find(MB(1, :) == MV);
    if size(TM, 2) ~= 1
        ii = 0;
        while (2*VAL.n) + 1 ~= ii
            ii  = ii + 1;
            if ~isnan(min(MSS.LFmin(ii, MB(3, TM)))) && length(TM) > 1
                TMM = TM(MSS.LFmin(ii, MB(3, TM)) == min(MSS.LFmin(ii, MB(3, TM))));
            end
        end
        if length(TMM) > 1
            pp = MB(3, TMM);
            ii = find(~isnan(MSS.LFmin(:, pp(1))), 1, 'last');
            TMM = TMM(find(MSS.LFmin(ii, MB(3, TMM)) ==...
                min(MSS.LFmin(ii, MB(3, TMM))), 1, 'last'));
        end
        [~, MI] = max(MB(3, TMM));
        MI = TMM(MI);
    end
    SMS(:, i) = MB(:, MI);
end
SMS   = flip(SMS, 2);
boxes = sort(SMS(3, :));
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function: Pareto
% Purpose:  Selection of potential optimal hyper-rectangles
%--------------------------------------------------------------------------
function POH = Pareto(VAL, MSS)
%--------------------------------------------------------------------------
% Calculate Euclidean Distatnces
aa = find(MSS.DD(1:VAL.I) ~= -1);
Euclid_dist = sum((VAL.Xmin(:, 1) - MSS.CC(:, aa)).^2, 1).^0.5;

% Identify potential optimal hyper-rectangles
S = aa(Find_poh_pareto(MSS.FE(aa), MSS.DD(aa), 1:length(aa), MSS, VAL, 1));
D = aa(Find_poh_pareto(Euclid_dist, MSS.DD(aa), 1:length(aa), MSS, VAL, 2));

% Find unique set of potential optimal hyper-rectangles
D(ismember(D, intersect(D, S))) = [];
D = sortrows([D; MSS.DD(D)].', 2).';
S = sortrows([S; MSS.DD(S)].', 2).';
POH = [S(1, :), D(1, :)];
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function: Find_poh_pareto
% Purpose: Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function boxes = Find_poh_pareto(fc, szes, indexs, MSS, VAL, phase)
%--------------------------------------------------------------------------
C   = unique(szes);
if C(1) == -1
   C(1) = []; 
end
aa  = find(MSS.DD(1:VAL.I) ~= -1);   
SMS = zeros(3, size(C, 2));
m_set = [fc; szes; indexs];
for i = 1: size(C, 2)                             % reduce m_set
    MB          = m_set(:, m_set(2, :) == C(i));
    [MV, MI]    = min(MB(1, :));
    TM          = find(MB(1, :) == MV);
    if size(TM, 2) ~= 1
        ii = 0;
        while (2*VAL.n) + 1 ~= ii
            ii  = ii + 1;
            if ~isnan(min(MSS.LFmin(ii, aa(MB(3, TM))))) && length(TM) > 1
                TM = TM(MSS.LFmin(ii, aa(MB(3, TM))) ==...
                    min(MSS.LFmin(ii, aa(MB(3, TM)))));
            else
                break;
            end
        end
        if length(TM) > 1
            pp = aa(MB(3, TM));
            ii = find(~isnan(MSS.LFmin(:, pp(1))), 1, 'last');
            TM = TM(find(MSS.LFmin(ii, aa(MB(3, TM))) ==...
                min(MSS.LFmin(ii, aa(MB(3, TM)))), 1, 'last'));
        end
        [~, MI]   = max(MB(3, TM));
        MI        = TM(MI);
    end
    SMS(:, i)   = MB(:, MI);
end
if phase == 3
    indexass = 1:length(C);
    for i = 1: size(C, 2)
        TM = indexass(SMS(1, indexass) == min(SMS(1, indexass)));
        PP = (find(SMS(1, indexass) == min(SMS(1, indexass))));
        if length(TM) > 1
            TMM = TM;
            ii = 0;
            while (2*VAL.n) + 1 ~= ii
                ii  = ii + 1;
                if ~isnan(min(MSS.LFmin(ii, aa(SMS(3, TMM))))) && length(TMM) > 1
                    TMM = TMM(MSS.LFmin(ii, aa(SMS(3, TMM))) ==...
                        min(MSS.LFmin(ii, aa(SMS(3, TMM)))));
                else
                    break;
                end
            end
            delet = setdiff(TM, TMM);
            delet = setdiff(delet, length(C));
            SMS(1, delet) = nan;
        end
        indexass(PP) = [];
    end
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
getopts(opts, 'maxits', 1000, 'maxevals', 100000, 'maxdeep', 1000,...
    'testflag', 0, 'tol', 0.01, 'showits', 1, 'dimension', 1,...
    'globalmin', 0, 'globalxmin', 0, 'ep', 1e-4, 'equation', '22d', 'poh', 'GL');

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
OPTI.ep        = ep;       % allowable relative error if f_reach is set
OPTI.equation  = equation; % equation option
OPTI.poh       = poh;      % POH selection scheme
VAL.equation   = OPTI.equation;
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function: Alocate
% Purpose: Create necessary data accessible to all workers
%--------------------------------------------------------------------------
function [VAL, MSS] = Alocate(OPTI, VAL)
%--------------------------------------------------------------------------
tic                                      % Mesure time
VAL.time = toc;                          % initial time
VAL.perror = 10;                        % initial perror

% alociate MAIN sets
MSS = struct('FF', zeros(1, OPTI.MAXevals+10000),...
             'FE', zeros(1, OPTI.MAXevals+10000),...
             'DD', zeros(1, OPTI.MAXevals+10000),...
             'LL', ones(VAL.n, OPTI.MAXevals+10000),...
             'CC', zeros(VAL.n, OPTI.MAXevals+10000),...
             'LFA', nan((2*VAL.n) + 1, OPTI.MAXevals+10000),...
             'LFmin', nan((2*VAL.n) + 1, OPTI.MAXevals+10000));

                    
if OPTI.G_nargout == 3
    VAL.history = zeros(OPTI.MAXits, 4); % allocating history
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function: Inits
% Purpose: Initialization of the HALRECT
%--------------------------------------------------------------------------
function [OPTI, VAL, MSS] = Inits(VAL, OPTI, Problem, MSS)
%--------------------------------------------------------------------------
VAL.itctr = 1;
[VAL.fminindex, VAL.I] = deal(1);
MSS.CC(:, 1) = ones(VAL.n, 1)/2;                  % initial midpoint
[MSS.FF(1), VAL.Fmin, MSS.FE(1)] = deal(feval(Problem.f, (abs(VAL.b -...
    VAL.a).*(MSS.CC(:, 1)) + VAL.a)));
MSS.LFA((2*VAL.n) + 1, 1) = 1;
VAL.Xmin(:, 1) = MSS.CC(:, 1);                        % initial point
MSS.DD(1) = norm(MSS.LL(:, 1),  2);
%--------------------------------------------------------------------------

% Check stop condition if global minima is known
if OPTI.TESTflag  == 1
    if OPTI.globalMIN ~= 0
        VAL.perror = 100*(VAL.Fmin - OPTI.globalMIN)/abs(OPTI.globalMIN);
    else
        VAL.perror = 100*VAL.Fmin;
    end
else
    VAL.perror = 2;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function  : Subdivision
% Purpose   : Calculate; new points, function values, lengths and indexes
%--------------------------------------------------------------------------
function [MSS, VAL] = Subdivisionas(VAL, Problem, MSS, POHa)
%-------------------------------------------------------------------------- 
% Find quantile of the diamteters set
for g = 1:size(POHa, 2)
    lsas = find(MSS.LL(:, POHa(1, g)) == max(MSS.LL(:, POHa(1, g))));
    pot = VAL.I;
    i = 0;
    partition_dim = 0;
    DIVIDE = POHa(1, g);
    stop = 3;
    while i ~= stop
        i = i + 1;
        if i ~= 1 && pot + i <= VAL.I
            indd = MSS.LFA(find(~isnan(MSS.LFA(:, pot + i))), pot + i);
            indf = MSS.LFA(find(~isnan(MSS.LFA(:, pot + i + 1))), pot + i + 1);
            dada = [min(MSS.FF(indd)), min(MSS.FF(indf))];
            mins = min(dada);
            inds = find(dada == mins, 1, "last");
            i = i + inds - 1;
            if mins == VAL.Fmin && stop <= (length(lsas))*2 + 1
                DIVIDE = pot + i;
                partition_dim = 0;
                stop = stop + 2;
                i = stop - 2;
            else
                partition_dim = 1;
            end
        end

        if partition_dim == 0
            max_L = max(MSS.LL(:, DIVIDE));
            ls = find(MSS.LL(:, DIVIDE) == max_L);
            if size(ls, 1) ~= 1
                DIR = max(abs(MSS.CC(:, VAL.fminindex) - MSS.CC(:, DIVIDE)), [], 2);
                ls = ls(find(DIR(ls) == max(DIR(ls)), 1, 'first'));
            end

            % Calculate side lengths
            MSS.LL(:, VAL.I + 1:VAL.I + 2) = MSS.LL(:, DIVIDE)*ones(1, 2);
            MSS.LL(ls, [VAL.I + 1, VAL.I + 2]) = max_L/2;

            % Calculate new points 

            MSS.CC(:, VAL.I + 1:VAL.I + 2) = MSS.CC(:, DIVIDE)*ones(1, 2);
            MSS.CC(ls, VAL.I + 1) = MSS.CC(ls, VAL.I + 1) - max_L/4;
            MSS.CC(ls, VAL.I + 2) = MSS.CC(ls, VAL.I + 2) + max_L/4;

            % Evaluate the objective function
            MSS.FF(VAL.I + 1) = feval(Problem.f, abs(VAL.b - VAL.a).*MSS.CC(:, VAL.I + 1) + VAL.a);
            MSS.FF(VAL.I + 2) = feval(Problem.f, abs(VAL.b - VAL.a).*MSS.CC(:, VAL.I + 2) + VAL.a);
            
            % Calculate diameters
            MSS.DD(VAL.I + 1:VAL.I + 2) = norm(MSS.LL(:, VAL.I + 1),  2);
            MSS.DD(DIVIDE) = -1;

            % Find index of the hyper-rectangles
            indd = MSS.LFA(find(~isnan(MSS.LFA(:, DIVIDE))), DIVIDE);

            jjj = [indd(MSS.CC(ls, indd) <= MSS.CC(ls, DIVIDE))];
            MSS.LFA(1:length(jjj), VAL.I + 1) = jjj;
            MSS.LFA((2*VAL.n) + 1, VAL.I + 1) = VAL.I + 1;
            MSS.LFmin(1:(length(jjj) + 1), VAL.I + 1) = sort(MSS.FF([jjj; VAL.I + 1]))';

            iii = [indd(MSS.CC(ls, indd) >= MSS.CC(ls, DIVIDE))];
            MSS.LFA(1:length(iii), VAL.I + 2) = iii;
            MSS.LFA((2*VAL.n) + 1, VAL.I + 2) = VAL.I + 2;
            MSS.LFmin(1:(length(iii) + 1), VAL.I + 2) = sort(MSS.FF([iii; VAL.I + 2]))';      
            
            BD = min([MSS.FF([jjj; VAL.I + 1]), MSS.FF(VAL.I + 1)]);
            BF = max([BD, MSS.FF(VAL.I + 1)]);
            if BD <= VAL.Fmin, VAL.Fmin = BD; end
          
            switch VAL.equation
                case '22a'
                    MSS.FE(VAL.I + 1) = MSS.FF(VAL.I + 1);
                case '22b'
                    MSS.FE(VAL.I + 1) = BD;
                case '22c'
                    MSS.FE(VAL.I + 1) = mean([MSS.FF([jjj; VAL.I + 1]), MSS.FF(VAL.I + 1)]);
                otherwise 
                    MSS.FE(VAL.I + 1) = 0.5*BD + 0.5*BF;
            end

            BP = min([MSS.FF([iii; VAL.I + 2]), MSS.FF(VAL.I + 2)]);
            BH = max([BP, MSS.FF(VAL.I + 2)]);
            if BP <= VAL.Fmin, VAL.Fmin = BP; end

            switch VAL.equation
                case '22a'
                    MSS.FE(VAL.I + 2) = MSS.FF(VAL.I + 2);
                case '22b'
                    MSS.FE(VAL.I + 2) = BP;
                case '22c'
                    MSS.FE(VAL.I + 2) = mean([MSS.FF([iii; VAL.I + 2]), MSS.FF(VAL.I + 2)]);
                otherwise
                    MSS.FE(VAL.I + 2) = 0.5*BP + 0.5*BH;
            end
            VAL.I = VAL.I + 2;
        end
    end
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function: Are we done
% Purpose:  Update minima value and check stopoing conditions
function [VAL, MSS] = Arewedone(OPTI, VAL, MSS)
%--------------------------------------------------------------------------
[VAL.Fmin, ~]  = min(MSS.FF(1:VAL.I));
VAL.fminindex = find(MSS.FF(1:VAL.I) == VAL.Fmin, 1,"last");
VAL.Xmin      = MSS.CC(:, VAL.fminindex);

if OPTI.showITS == 1                % Show iteration stats
    VAL.time = toc;
    fprintf(...
    'Iter: %4i   f_min: %15.10f    time(s): %10.05f    fn evals: %8i\n',...
        VAL.itctr - 1, VAL.Fmin, VAL.time, VAL.I);
end

if OPTI.TESTflag == 1               % Check for stop condition
    if OPTI.globalMIN ~= 0          % Calculate error if globalmin known
        VAL.perror = 100*(VAL.Fmin - OPTI.globalMIN)/abs(OPTI.globalMIN);
    else
        VAL.perror = 100*VAL.Fmin;
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

if VAL.I > OPTI.MAXevals            % Have we exceeded the maxevals?
    disp('Exceeded max fcn evals. Increase maxevals');
    VAL.perror = -10;
end

if max(max(MSS.LL)) + 1 > OPTI.MAXdeep  % Have we exceeded max deep?
    disp('Exceeded Max depth. Increse maxdeep');
    VAL.perror = -10;
end

if OPTI.G_nargout == 3              % Store History
    VAL.history(VAL.itctr,1) = VAL.itctr;
    VAL.history(VAL.itctr,2) = VAL.I;
    VAL.history(VAL.itctr,3) = VAL.Fmin;
    VAL.history(VAL.itctr,4) = VAL.time;
end
% Update iteration number
if VAL.perror > OPTI.TOL
    VAL.itctr = VAL.itctr + 1;
end
%--------------------------------------------------------------------------
return

%--------------------------------------------------------------------------
% Function:  Find_po
% Purpose:  Return list of PO hyperrectangles
%--------------------------------------------------------------------------
function S = Find_po(fc, f_min, epsilon, szes)
% Find all rects on hub
E          = max(epsilon*abs(f_min), 1E-8);
[~, i_min] = min((fc - f_min + E)./szes);
d          = unique(szes);
idx        = find(d == szes(i_min));
jj         = length(d);
d_min      = zeros(1, length(d));

for i = 1:jj
    d_min(i) = min(fc(szes == d(i)));
end

Su = cell(1, jj);
for i = idx:jj
    idx2 = find((fc == d_min(i)) & (szes == d(i)), 1, 'first');
    Su{i} = idx2;
end
S_1 = cell2mat(Su);
Su = cell(1, jj);
if length(d) - idx > 1
    a1 = szes(i_min);
    b1 = fc(i_min);
    a2 = d(length(d));
    b2 = d_min(length(d));
    % The line is defined by: y = slope*x + const
    slope = (b2-b1)/(a2-a1);
    const = b1 - slope*a1;
    for i = 1 : length(S_1)
        j = S_1(i); 
        if fc(j) <= slope*szes(j) + const + 1E-6
            Su{i} = j;
        end
    end
    S_2 = cell2mat(Su);
    if isempty(S_2)
        S_2 = S_1;
    end
    % S_2 now contains all points in S_1 which lies on or below the line
    % Find the points on the convex hull defined by the points in S_2
    xx = szes(S_2);
    yy = fc(S_2);
    h = conhull(xx, yy); % conhull is an internal subfunction
    S_3 = S_2(h);
else
    S_3 = S_1;
end
S = S_3;
return

%--------------------------------------------------------------------------
% Function:  conhull
% Purpose:  conhull returns all points on the convex hull.
%--------------------------------------------------------------------------
function h = conhull(x, y)
% conhull returns all points on the convex hull, even redundant ones.
format long;
x = x(:);
y = y(:);
xyAR = [x y];
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