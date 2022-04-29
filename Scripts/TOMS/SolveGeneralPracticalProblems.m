% -------------------------------------------------------------------------
% Script: SolveGeneralPracticalProblems
%
% Author 1: Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2: Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
%
% Created on: 04/20/2022
%
% Purpose 
% This code is prepared for minimization of six practical problems with
% general constraints.
%
% Source 
% Stripinis, L., Paulavicius, R.: DIRECTGO: A new DIRECT-type MATLAB toolbox
% for derivative-free global optimization (2022)
%
% OUTPUTS
% SolveGeneralPracticalProblems_results.mat
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
parts{end + 1} = 'Engineering'; 
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
% parts{end + 1} = 'Engineering'; 
% parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 

% Load path:
% -------------------------------------------------------------------------
parts = strsplit(pwd, filesep); parts{end} = 'Algorithms';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path);  
parts = strsplit(pwd, filesep); parts = parts(1:end-1);
parts{end} = 'Results'; parts{end + 1} = 'TOMS';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
load('Results_on_constrained_engineering_problems_(eps = 0.01).mat');

% Results:
% -------------------------------------------------------------------------
SolveGeneralPracticalProblems_results = Results4(2:size(Results4, 1), 1:3);
SolveGeneralPracticalProblems_results{1, 4} = "Iterations";
SolveGeneralPracticalProblems_results{1, 5} = "Evaluations";
SolveGeneralPracticalProblems_results{1, 6} = "Seconds";

% Load options:
% -------------------------------------------------------------------------
opts.maxevals = 10000;      % Maximum number of function evaluations.
opts.maxits   = 10000;      % Maximum number of iterations.
opts.maxdeep  = 10000;      % Maximum number of side divisions.
opts.testflag = 1;          % Testing algorithms with known globalmin vcalue.
opts.showits  = 1;          % Print iteration stats.
opts.tol      = 0.01;       % Allowable relative error if globalmin is set.
opts.ept      = 0.00000001; % Tollerance for constrains

% Load options:
% -------------------------------------------------------------------------
for h = 2:size(SolveGeneralPracticalProblems_results, 1)
    % Number of variables
    opts.dimension = SolveGeneralPracticalProblems_results{h, 3}; 

    % Select problem
    Problem.f = SolveGeneralPracticalProblems_results{h, 2};

    % Solve problem
    [Fmin, Xmin, history] = dDirect_GLc(Problem, opts);

    % Record algorithm performance results  
    SolveGeneralPracticalProblems_results{h, 4} =  history(end, 1);
    SolveGeneralPracticalProblems_results{h, 5} =  history(end, 2);
    SolveGeneralPracticalProblems_results{h, 6} =  history(end, 4);
end

% Store results:
% -------------------------------------------------------------------------
filename = 'SolveGeneralPracticalProblems_results.mat'; save(filename);