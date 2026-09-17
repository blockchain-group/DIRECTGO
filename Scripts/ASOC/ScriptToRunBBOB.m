% -------------------------------------------------------------------------
% Script: ScriptToRunBBOB.m
% Authors: Linas Stripinis, Remigijus Paulavicius
% Created: 09/09/2026
%
% Purpose:
% Benchmark dX_DTC_GL on the 24 noiseless BBOB functions over multiple
% dimensions and instances using IOH/IOHexperimenter.
%
% Output:
% Results.mat
%
% Requirements:
% - MATLAB with Python support 
% - Python with the IOH package installed (pip install ioh) 
%
% References:
% [1] N. Hansen, S. Finck, R. Ros, and A. Auger,
%     "Real-Parameter Black-Box Optimization Benchmarking 2009:
%     Noiseless Functions Definitions," INRIA, Tech. Rep. RR-6829, 2009.
%
% [2] IOHprofiler, "IOHexperimenter: Experimenter for Iterative
%     Optimization Heuristics," https://github.com/IOHprofiler/IOHexperimenter
% -------------------------------------------------------------------------

clear; clc;
addpath('Algorithms');

% Check Python configuration
pe = pyenv;
if pe.Status == "NotLoaded" && strlength(pe.Executable) == 0
    error('Python is not configured. Configure it using pyenv.');
end

% Check IOH package 
try 
    ioh = py.importlib.import_module('ioh'); 
catch
    error('Python package "ioh" is required. Install it in the Python environment used by MATLAB: pip install ioh'); 
end

%% Experimental setup
MaxEvals  = 1e5;               % Evaluations per dimension
MaxIts    = 1e6;               % Maximum iterations
Dimensions = [2, 3, 5, 10, 20];
Instances  = 5;
Error      = 1e-4;

%% Results
DIRECTGOLib_Results = {"Nr.", "Problem name", "Instance", "Dimension", "Fmin", "Xmin", "History", "Fbest", "Xbest"};

%% Run experiments
for h = 1:24
    for i = 1:Instances
        for n = Dimensions

            global history ii mma xx %#ok<*TLEV,*GVMIS>

            % Algorithm options
            opts.dimension = n;
            opts.maxevals  = MaxEvals*n;
            opts.maxits    = MaxIts;
            opts.testflag  = 2;

            % Bounds
            Bounds = [-5*ones(n,1), 5*ones(n,1)];

            % BBOB problem
            fhd = ioh.get_problem(py.int(h),py.int(i),py.int(n),ioh.ProblemClass.BBOB);

            f_opt        = double(fhd.optimum.y);
            x_opt        = double(fhd.optimum.x.tolist())';
            opts.ftarget = f_opt + Error;
            Problem.f    = @(x) IOHfun(x, fhd);

            % Initialize evaluation history
            history = [0, 0, inf, 0];
            ii = 0;
            xx = [];
            mma = opts.maxevals;

            % Run algorithm
            tic;
            [fbest, xbest, ~] = dX_DTC_GL(Problem, opts, Bounds);

            % Store results
            k = size(DIRECTGOLib_Results, 1) + 1;
            DIRECTGOLib_Results(k,:) = {h, string(fhd.meta_data.name), i, n, f_opt, x_opt, history, fbest, xbest};

            fprintf('Function: %2i | Dimension: %2i | Instance: %2i\n', h, n, i);
        end
    end
end

%% Save results
save('Results.mat', 'DIRECTGOLib_Results');

function y = IOHfun(x, fhd)
    global history ii mma xx
    
    y = cellfun(@double, cell(fhd(py.list({x}))));
    ii = ii + 1;
    
    if y < history(end,3) && ii <= mma
        xx = x;
        history(end+1,:) = [size(history,1)-1, ii, y, toc];
    elseif ii == mma
        history(end+1,:) = [size(history,1)-1, ii, history(end,3), toc];
    end
end