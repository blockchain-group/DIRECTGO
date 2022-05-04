% -------------------------------------------------------------------------
% Script: SolveHALRECT
%
% Author 1: Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2: Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
%
% Created on: 04/27/2022
%
% Purpose 
% This code is prepared for minimization of 81 test problems with
% box constraints.
%
% Source 
% Stripinis, L., Paulavicius, R.: DIRECTGO: A new DIRECT-type MATLAB toolbox
% for derivative-free global optimization (2022)
%
% OUTPUTS
% HALRECT_DIRECTGOLib_results.mat
%--------------------------------------------------------------------------

% Download DIRECTGOLib v1.1
% Stripinis, L., Paulavicius, R.: DIRECTGOLib - DIRECT Global Optimization 
% test problems Library. Zenodo (2022). doi.org/10.5281/zenodo.6491863
% -------------------------------------------------------------------------
if not(isfolder('DIRECTGOLib-1.1'))
    fullURL = 'https://github.com/blockchain-group/DIRECTGOLib/archive/refs/tags/v1.1.zip';
    filename = 'DIRECTGOLib.zip';
    websave(filename, fullURL);
    unzip('DIRECTGOLib.zip');
end
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-1.1';
parts{end + 1} = 'Box';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path);

% Download newest version of DIRECTGOLib
% Stripinis, L., Paulavicius, R.: DIRECTGOLib - DIRECT Global Optimization 
% test problems Library. Zenodo (2022). doi.org/10.5281/zenodo.6491951
% -------------------------------------------------------------------------
% if not(isfolder('DIRECTGOLib-main'))
%     fullURL = 'https://github.com/blockchain-group/DIRECTGOLib/archive/refs/heads/main.zip';
%     filename = 'DIRECTGOLib.zip';
%     websave(filename, fullURL);
%     unzip('DIRECTGOLib.zip'); 
% end
% parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-main'; 
% parts{end + 1} = 'Box'; 
% parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 

% Load path:
% -------------------------------------------------------------------------
parts = strsplit(pwd, filesep); parts = parts(1:end-1); 
parts{end} = 'Algorithms'; parent_path = strjoin(parts(1:end), filesep);
addpath(parent_path); 
parts = strsplit(pwd, filesep); parts = parts(1:end-1);
parts{end} = 'Results'; parts{end + 1} = 'COA';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
load('HALRECT_results.mat');
% -------------------------------------------------------------------------

% Load options:
% -------------------------------------------------------------------------
perturb       = 0;       % Percentage perturbation of bound constraints
% -------------------------------------------------------------------------
opts.maxevals = 10000;   % Maximum number of function evaluations.
opts.maxits   = 10000;   % Maximum number of iterations.
opts.maxdeep  = 10000;   % Maximum number of side divisions.
opts.testflag = 1;       % Testing algorithms with known globalmin vcalue.
opts.showits  = 1;       % Print iteration stats.
opts.tol      = 0.01;    % Allowable relative error if globalmin is set.
opts.ep       = 0.0001;  % global/local weight parameter.
opts.poh      = 'GL';    % Choose the the selection of POH scheme:
                         % 'GL' - Two-step-based Pareto selection (default) 
                         % 'IA' - Improved Aggressive selection strategy 
                         % 'IO' - Improved Original selection strategy
opts.equation = '22d';   % Selection is performed using one of the four 
                         % posible models: '22a', '22b', '22c', '22d'(default)

% Results:
% -------------------------------------------------------------------------
EmpiricalDIRECTGOLib = Results(1:size(Results, 1), 1:5);
EmpiricalDIRECTGOLib{1, 6} = "Evaluations";

% Solve problems:
% -------------------------------------------------------------------------
for h = 2:size(EmpiricalDIRECTGOLib, 1)
    % Number of variables
    opts.dimension = EmpiricalDIRECTGOLib{h, 3};   

    % Select problem
    Problem.f = EmpiricalDIRECTGOLib{h, 2};   

    % Select minima
    y = feval(Problem.f);
    opts.globalmin = y.fmin(opts.dimension);
    opts.globalXMIN = y.xmin(opts.dimension);

    % Select bound constraints
    x_lower = EmpiricalDIRECTGOLib{h, 4};
    x_upper = EmpiricalDIRECTGOLib{h, 5};
    shifts = abs(x_lower - x_upper)*perturb;
    x_lower = min([opts.globalXMIN, x_lower + shifts], [], 2);
    x_upper = x_upper + shifts;
    bounds = [x_lower, x_upper];
 
    % Solve problem
    [~, ~, history] = HALRECT(Problem, opts, bounds);

    % Record algorithm performance results   
    EmpiricalDIRECTGOLib{h, 6} =  history(end, 2);   
end

% Store results:
% -------------------------------------------------------------------------
filename = 'HALRECT_DIRECTGOLib_results.mat'; save(filename);