% -------------------------------------------------------------------------
% Script: SolveDIRECTGOlib
%
% Author 1: Linas Stripinis          (linas.stripinis@mif.vu.lt)
% Author 2: Remigijus Paualvicius    (remigijus.paulavicius@mif.vu.lt)
%
% Created on: 28/07/2023
%
% Purpose 
% This code is prepared for minimization of 324 test problems with
% box constraints.
%
% OUTPUTS
% Results.mat
%--------------------------------------------------------------------------
% Download newest version of DIRECTGOLib
%  Stripinis, L., Kudela, J., Paulavicius, R.: Directgolib - direct global 
%  optimization test problems library (2023). Pre-release v2.0, URL 
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
parts{end + 1} = 'box'; parts{end + 1} = 'ABS';  
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-main'; 
parts{end + 1} = 'box'; parts{end + 1} = 'BBOB';  
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-main'; 
parts{end + 1} = 'box'; parts{end + 1} = 'CEC';  
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-main'; 
parts{end + 1} = 'box'; parts{end + 1} = 'Classical';  
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
parts = strsplit(pwd, filesep); parts{end + 1} = 'DIRECTGOLib-main'; 
parts{end + 1} = 'box'; parts{end + 1} = 'Layeb';  
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
% -------------------------------------------------------------------------
% Load other path:
% -------------------------------------------------------------------------
parts = strsplit(pwd, filesep); parts = parts(1:end-1); 
parts{end} = 'Algorithms'; parent_path = strjoin(parts(1:end), filesep);
addpath(parent_path); 
parts = strsplit(pwd, filesep); parts = parts(1:end-1);
parts{end} = 'Results'; parts{end + 1} = 'MPC';
parent_path = strjoin(parts(1:end), filesep); addpath(parent_path); 
load('DIRECTGOLib_results.mat');
% -------------------------------------------------------------------------
%% Experimental Setup
% Maximum number of function evaluations
MaxEvals = 1e1;  % should also get multiplied by dimension later on

% Maximum number of iterations
MaxIts = 10^6; 

% Allowable relative error if globalmin is set
Perror = 1e-2;

%% Loop over all prepared test problems:
for h = 2:size(DIRECTGOLib_Results, 1) 
    % Extract info from the problem:
    [dim, fun, inst, xL, xU, Fmin, Xmin] = ExtractingInfo(DIRECTGOLib_Results, h);
    
    % An example of the use of the DIRECT algorithm:
    opts = struct(); Problem = struct();
    opts.dimension  = dim;
    opts.globalmin  = Fmin;
    opts.globalxmin = Xmin;
    opts.maxevals   = MaxEvals*dim;
    opts.maxits     = MaxIts; 
    opts.tol        = Perror;
    opts.showits    = 1;
    opts.testflag   = 1;
    opts.population = 100;
    Bounds          = [xL, xU];

    M = eye(dim);
    [shift_min, shift_max] = compute_shift_bounds(M,Xmin,xL,xU,inst);
    rng(inst,'twister'); %seed inst for shift
    shift = (shift_min + 0.1*(shift_max - shift_min).*rand(dim,1));

    getInfo = feval(fun);
    if getInfo.libraries(9) == 1 || getInfo.libraries(10) == 1
        func = str2func(['@(x, shift, M) ',fun,'(x, shift, M)']);
        fun_rot = @(x) func(x, shift, M);
    else
        func = str2func(['@(x) ',fun,'(x)']);
        temp_vec = -M*shift;
        fun_rot = @(x) func(min(max(M*x + temp_vec,xL),xU));
    end

    Problem.f = fun_rot; 
    [fbest, xatmin, history] = dDirect_GL(Problem, opts, Bounds);

    DIRECTGOLib_Results{h, 8} = history; %#ok<*SAGROW>
    DIRECTGOLib_Results{h, 9} = fbest;
    DIRECTGOLib_Results{h, 10} = xatmin;
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
    Xmin = getInfo.xmin(dim);
end

function [t_max] = shiftLP(M,Xmin,xL,xU,shift_dir)
    options = optimoptions('linprog','Display','none');
    Aeq = [M,-M*shift_dir]; beq = Xmin;
    temp_max = linprog([0*Xmin;-1],[],[],Aeq,beq,[xL;-inf], [xU;inf],options);
    t_max = temp_max(end);
end

function [shift_min, shift_max] = compute_shift_bounds(M,Xmin,xL,xU,inst)
    dim = length(Xmin);
    rng(inst,'twister');
    shift_dir = randn(dim,1);
    shift_min = (0*shift_dir);
    rng(inst,'twister');
    [t_max] = shiftLP(M,Xmin,xL,xU,shift_dir);
    shift_max = t_max*shift_dir;
    rng(inst,'twister');
end