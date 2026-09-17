% -------------------------------------------------------------------------
% Script: ScriptToRunDIRECTGOLib.m
% Authors: Linas Stripinis, Remigijus Paulavicius
% Created: 09/09/2026
%
% Purpose:
% Run the selected DIRECTGOLib v2.0 box-constrained test instances using
% the chosen DIRECT-type algorithm and store the performance results.
%
% Output:
% Results.mat
%
% References:
% [1] Stripinis, L., Kudela, J., Paulavicius, R.: Directgolib - direct global 
%     optimization test problems library (2024). Pre-release v2.0, URL 
%     https://github.com/blockchain-group/DIRECTGOLib.
%
% [2] L. Stripinis, J. Kůdela, R. Paulavičius, Two novel instance selection
%     methods combining algorithm performance and landscape analysis: A
%     comparative study in continuous optimization, IEEE Transactions on
%     Cybernetics 56 (3) (2026) 1202–1215. doi:10.1109/TCYB.2025.3625095.
% -------------------------------------------------------------------------

%% Experimental setup
MaxEvals = 1e5;        % Evaluations per dimension
MaxIts   = 1e6;        % Maximum iterations
Error    = 1e-4;       % error
BS       = 1;          % Selected subset of problems using ISMs [2]

%% Paths
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
addpath('Algorithms','DIRECTGOLib-main');

%% Load test instances
switch BS
    case 1, load('DIRECTGOLib_settings_cs.mat');
    case 2, load('DIRECTGOLib_settings_md.mat');
    case 3, load('DIRECTGOLib_settings_rb.mat');
    case 4, load('DIRECTGOLib_settings_ped.mat');
end

%% Run experiments
for h = 2:size(DIRECTGOLib_Results,1)

        clear functions
        global history ii mma xx %#ok<*TLEV,*GVMIS>

        % Problem data
        [dim, fun, xL, xU, Fmin, M, shift] = ExtractingInfo(DIRECTGOLib_Results,h);

        fhd = compute_function(xL,xU,fun,M,shift);

        % Algorithm options
        opts.dimension = dim;
        opts.maxevals  = MaxEvals*dim;
        opts.maxits    = MaxIts;
        opts.showits   = 1;
        opts.testflag  = 2;
        opts.ftarget   = Fmin + 1e-4;

        Bounds    = [xL,xU];
        Problem.f = @(x) funkcija(fhd,x);

        % Initialize evaluation history
        history = [0,0,inf,0];
        ii = 0;
        xx = [];
        mma = opts.maxevals;

        % Run algorithm
        tic
        [fbest,xbest,~] = dX_DTC_GL(Problem,opts,Bounds);

        % Store results
        DIRECTGOLib_Results{h,7}  = Fmin; %#ok<*SAGROW>
        DIRECTGOLib_Results{h,8}  = history;
        DIRECTGOLib_Results{h,9}  = fbest;
        DIRECTGOLib_Results{h,10} = xbest;
        DIRECTGOLib_Results{h,11} = history(end,2);
end
save('Results.mat', 'DIRECTGOLib_Results');

function [dim,fun,xL,xU,Fmin,M,shift] = ExtractingInfo(R,h)
    
    dim = R{h,3};
    fun = R{h,2};
    
    info = feval(fun);
    
    xL   = info.xl(dim);
    xU   = info.xu(dim);
    Fmin = info.fmin(dim);
    
    M     = R{h,5};
    shift = R{h,6};
end


function fhd = compute_function(xL,xU,fun,M,shift)
    
    info = feval(fun);
    xM   = (xL+xU)/2;
    
    if info.libraries(9) ~= 1 && info.libraries(10) ~= 1
    
        f = str2func(['@(x)',fun,'(x)']);
        t = -M*shift - M*xM + xM;
    
        fhd = @(x) f(min(max(M*x+t,xL),xU));
    
    else
    
        f = str2func(['@(x,shift,M)',fun,'(x,shift,M)']);
        fhd = @(x) f(x,shift,M);
    
    end
end


function y = funkcija(fhd,x)
    
    global history ii mma xx
    
    y  = min(feval(fhd,x),1e300);
    ii = ii + 1;
    
    if y < history(end,3) && ii <= mma
        xx = x;
        history(end+1,:) = [size(history,1)-1, ii, y, toc];
    
    elseif ii == mma
        history(end+1,:) = [size(history,1)-1, ii, history(end,3), toc];
    end
end