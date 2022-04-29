% -------------------------------------------------------------------------
% Script: SolveBoxProblems
%
% Author 1: Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2: Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
%
% Created on: 04/26/2022
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
% SolveBoxProblems_results.mat
%--------------------------------------------------------------------------

% Download DIRECTGOLib v1.0
% Stripinis, L., Paulavicius, R.: DIRECTGOLib - DIRECT Global Optimization 
% test problems Library. Zenodo (2022). doi.org/10.5281/zenodo.6491863
% -------------------------------------------------------------------------
if not(isfolder('DIRECTGOLib-1.0'))
    fullURL = 'https://github.com/blockchain-group/DIRECTGOLib/archive/refs/tags/v1.0.zip';
    filename = 'DIRECTGOLib.zip';
    websave(filename, fullURL);
    unzip('DIRECTGOLib.zip');
end
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-1.0';
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
parts{end} = 'Results'; parts{end + 1} = 'TOMS';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
load('Results_on_box_constrained_problems_(eps = 0.01).mat');
% -------------------------------------------------------------------------

% Load options:
% -------------------------------------------------------------------------
opts.maxevals = 10000;   % Maximum number of function evaluations.
opts.maxits   = 10000;   % Maximum number of iterations.
opts.maxdeep  = 10000;   % Maximum number of side divisions.
opts.testflag = 1;       % Testing algorithms with known globalmin vcalue.
opts.showits  = 1;       % Print iteration stats.
opts.tol      = 0.01;    % Allowable relative error if globalmin is set.

% Results:
% -------------------------------------------------------------------------
SolveBoxProblems_results = Results1(2:size(Results1, 1), 1:5);
SolveBoxProblems_results{1, 6} = "Iterations";
SolveBoxProblems_results{1, 7} = "Evaluations";
SolveBoxProblems_results{1, 8} = "Seconds";

% Solve problems:
% -------------------------------------------------------------------------
for h = 2:size(SolveBoxProblems_results, 1)
    % Number of variables
    opts.dimension = SolveBoxProblems_results{h, 3};   

    % Select problem
    Problem.f = SolveBoxProblems_results{h, 2};   

    % Select minima
    y = feval(Problem.f);
    opts.globalmin = y.fmin(opts.dimension);

    % Select bound constraints
    x_lower = SolveBoxProblems_results{h, 4};
    x_upper = SolveBoxProblems_results{h, 5};
    bounds = [x_lower, x_upper];
 
    % Solve problem
    [~, ~, history] = dDirect_GL(Problem, opts, bounds);

    % Record algorithm performance results  
    SolveBoxProblems_results{h, 6} =  history(end, 1);   
    SolveBoxProblems_results{h, 7} =  history(end, 2);   
    SolveBoxProblems_results{h, 8} =  history(end, 4);
end

% Store results:
% -------------------------------------------------------------------------
filename = 'SolveBoxProblems_results.mat'; save(filename);