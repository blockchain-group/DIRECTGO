% -------------------------------------------------------------------------
% Script: ScriptToRun
%
% Author 1: Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2: Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
%
% Created on: 01/09/2024
%
% Purpose 
% Loop through 24 BBOB test functions with different shifts and rotations
%
% OUTPUTS
% Results.mat
%--------------------------------------------------------------------------
% Download newest version of DIRECTGOLib
%  Stripinis, L., Kudela, J., Paulavicius, R.: Directgolib - direct global 
%  optimization test problems library (2024). Pre-release v2.0, URL 
%  https://github.com/blockchain-group/DIRECTGOLib.
% -------------------------------------------------------------------------
clear;clc;
if not(isfolder('DIRECTGOLib-main'))
    fullURL = 'https://github.com/blockchain-group/DIRECTGOLib/archive/refs/heads/main.zip';
    filename = 'DIRECTGOLib.zip';
    websave(filename, fullURL);
    unzip('DIRECTGOLib.zip'); 
end
% Load path:
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-main'; 
parts{end + 1} = 'box'; parts{end + 1} = 'BBOB';  
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
addpath('GENDIRECT')

PR = dir(parent_path);
names = {PR(3:26).name};
names_truncated = cellfun(@(x) x(1:end-2), names, 'UniformOutput', false);
% -------------------------------------------------------------------------
%% Experimental Setup
% Maximum number of function evaluations
MaxEvals = 1e5;  

% Maximum time limit
MaxTime = 1e10; 

% Maximum number of iterations
MaxIts = 10^6; 

% Allowable relative error if globalmin is set
Perror = 1e-8;

% Considered dimensions for scalable test functions
Dimensions = [2, 3, 5, 10, 20];

% Considered number of instances - first one without shift/rotation
Instances = 10;

%% Test instances:  
DIRECTGOLib_Results{1, 1} = "Nr."; DIRECTGOLib_Results{1, 2} = "Problem name";
DIRECTGOLib_Results{1, 3} = "Dimension"; DIRECTGOLib_Results{1, 4} = "Instance";
DIRECTGOLib_Results{1, 5} = "Bounds"; DIRECTGOLib_Results{1, 6} = "Xmin"; 
DIRECTGOLib_Results{1, 7} = "Fmin"; DIRECTGOLib_Results{1, 8} = "History"; 
DIRECTGOLib_Results{1, 9} = "Fbest"; DIRECTGOLib_Results{1, 10} = "Xbest";

% Prepare test functions:
for h = 1:length(names_truncated)
    for j = 1:length(Dimensions)
        ii = size(DIRECTGOLib_Results, 1);
        for jj=1:Instances
            DIRECTGOLib_Results{ii + jj, 1} = ii;
            DIRECTGOLib_Results{ii + jj, 2} = names_truncated{h};
            DIRECTGOLib_Results{ii + jj, 3} = Dimensions(j);
            DIRECTGOLib_Results{ii + jj, 4} = jj;
        end
    end
end

%% Loop over all prepared test problems:
for h = 2:size(DIRECTGOLib_Results, 1) 
    clear functions %#ok<*CLFUNC>
    
    % Extract info from the problem:
    [dim, fun, inst, xL, xU, Fmin, Xmin] = ExtractingInfo(DIRECTGOLib_Results, h);
    func = str2func(['@(x, inst) ', fun, '(x, inst)']);
    fun = @(x) func(x, inst);
    f_goal = Fmin + 1e-8;

    % Solve the test problem
    alg = GENDIRECT();
    alg.Problem.f = fun;
    alg.Problem.n = dim;
    alg.Problem.info = false;
    alg.Problem.xL = xL;
    alg.Problem.xU = xU;
    alg.Problem.fgoal = f_goal;

    alg.Partitioning.Strategy = 'DTC';
    alg.Partitioning.SubSides = 'All';
    alg.Partitioning.AggrFuncVal = 'Midpoint';

    alg.Selection.Strategy = 'Original';
    alg.Selection.CandMeasure = 'Diagonal';
    alg.Selection.EqualCand = 'All';
    alg.Selection.SolRefin = 'Off';
    alg.Selection.ControlEp = 'Off';
    alg.Selection.Ep = 1e-4;
    alg.Selection.GloballyBiased = 'Off';
    alg.Selection.TwoPhase = 'On';

    alg.Hybridization.Strategy = 'Off';
    alg.Hybridization.MaxEvaluations = 1000*dim;
    alg.Hybridization.MaxIterations = 1000*dim;
    alg.Hybridization.LocalSearch = 'sqp';

    alg.optParam.maxevals = MaxEvals*dim;
    alg.optParam.maxits = MaxEvals*dim;
    alg.optParam.showits = 1;

    alg = alg.solve;

    % Record algorithm performance results
    DIRECTGOLib_Results{h, 5} = [xL, xU];
    DIRECTGOLib_Results{h, 6} = Xmin;
    DIRECTGOLib_Results{h, 7} = Fmin;
    DIRECTGOLib_Results{h, 8} = alg.History;
    DIRECTGOLib_Results{h, 9} = alg.Fmin;
    DIRECTGOLib_Results{h, 10} = alg.Xmin;
end													

%% Store results:
save('Results','DIRECTGOLib_Results')

%% Extract info from the problem:
function [dim, fun, inst, xL, xU, Fmin, Xmin] = ExtractingInfo(DIRECTGOLib_Results, h)
    % Select dimension of the test problem
    dim = DIRECTGOLib_Results{h, 3}; 

    % Select the test problem
    fun = DIRECTGOLib_Results{h, 2}; 
    
    % Select instance
    inst = DIRECTGOLib_Results{h, 4}; 

    % Extract information from the test problem
    getInfo = feval(fun);

    % Bound constraints
    xL = getInfo.xl(dim);
    xU = getInfo.xu(dim);

    % Solution
    Fmin = getInfo.fmin(dim);
    Xmin = getInfo.xmin(dim, inst);
end