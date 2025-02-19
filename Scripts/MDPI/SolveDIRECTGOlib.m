% -------------------------------------------------------------------------
% Script: SolveDIRECTGOlib
%
% Author 1: Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2: Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
%
% Created on: 06/1/2023
%
% Purpose 
% This code is prepared for minimization of 67 test problems with
% linear constraints.
%
% Source 
% Stripinis, L., Paulavicius, R.: DIRECTGO: Novel Algorithm for Linearly 
% Constrained Derivative Free Global Optimization of Lipschitz Functions (2023)
%
% OUTPUTS
% EmpiricalDIRECTGOLib_results.mat
%--------------------------------------------------------------------------

% Download DIRECTGOLib v1.3
% Stripinis, L., Paulavicius, R.: DIRECTGOLib - DIRECT Global Optimization 
% test problems Library. Zenodo (2022). doi.org/10.5281/zenodo.6491863
% -------------------------------------------------------------------------
if not(isfolder('DIRECTGOLib-1.3'))
    fullURL = 'https://github.com/blockchain-group/DIRECTGOLib/archive/refs/tags/v1.3.zip';
    filename = 'DIRECTGOLib.zip';
    websave(filename, fullURL);
    unzip('DIRECTGOLib.zip');
end
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-1.3';
parts{end + 1} = 'Linear';
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
% parts{end + 1} = 'Linear'; 
% parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 

% Load path:
% -------------------------------------------------------------------------
parts = strsplit(pwd, filesep); parts = parts(1:end-1); 
parts{end} = 'Algorithms'; parent_path = strjoin(parts(1:end), filesep);
addpath(parent_path); 
parts = strsplit(pwd, filesep); parts = parts(1:end-1);
parts{end} = 'Results'; parts{end + 1} = 'MDPI';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
load('DIRECTGOLib_results.mat');
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
EmpiricalDIRECTGOLib = Results(1:size(Results, 1), 1:5);
EmpiricalDIRECTGOLib{1, 6} = "Evaluations";
EmpiricalDIRECTGOLib{1, 7} = "Minima";

% Solve problems:
% -------------------------------------------------------------------------
for h = 65:size(EmpiricalDIRECTGOLib, 1)
    % Number of variables
    opts.dimension = EmpiricalDIRECTGOLib{h, 3};   

    % Select problem
    Problem.f = EmpiricalDIRECTGOLib{h, 2};   
     
    % Select minima
    y = feval(Problem.f);
    opts.globalmin = y.fmin(opts.dimension);
    
    % Select constraint functions
    Problem.constraint = y.confun; 
    
    % Select bound constraints
    x_lower = EmpiricalDIRECTGOLib{h, 4};
    x_upper = EmpiricalDIRECTGOLib{h, 5};
    bounds = [x_lower, x_upper];
 
    % Solve problem
    [~, ~, history] = mBIRECTv_GL(Problem, opts, bounds);

    % Record algorithm performance results   
    EmpiricalDIRECTGOLib{h, 6} =  history(end, 2);   
    EmpiricalDIRECTGOLib{h, 7} =  history(end, 3); 
end

% Store results:
% -------------------------------------------------------------------------
filename = 'EmpiricalDIRECTGOLib_results.mat'; save(filename);