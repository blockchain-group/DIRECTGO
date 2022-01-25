% Load path:
% -------------------------------------------------------------------------
parts = strsplit(pwd, filesep);
parts{end} = 'Algorithms';
parent_path = strjoin(parts(1:end), filesep);
addpath(parent_path); 
parts{end} = 'DIRECTlib';
parts{end + 1} = 'Box';
parent_path = strjoin(parts(1:end), filesep);
addpath(parent_path);  
% -------------------------------------------------------------------------

% Load options:
% -------------------------------------------------------------------------
opts.maxevals = 10000; % Maximum number of function evaluations.
opts.maxits   = 10000; % Maximum number of iterations.
opts.maxdeep  = 10000;   % Maximum number of side divisions.
opts.testflag = 1;       % Testing algorithms with known globalmin vcalue.
opts.showits  = 1;       % Print iteration stats.
opts.tol      = 0.01;    % Allowable relative error if globalmin is set.
% -------------------------------------------------------------------------

% Test problems where the number of variables can be any integer n:
opts.dimension = 2;      % Specify number of variables
% -------------------------------------------------------------------------
BoxProblems ={'Ackley', 'Alpine', 'Csendes', 'Dixon_and_Price',...
    'Griewank', 'Levy', 'Michalewicz', 'Perm', 'Permdb', 'Qing',...
    'Rastrigin', 'Rosenbrock', 'Rotated_Hyper_Ellipsoid', 'Schwefel',...
    'Sphere', 'Styblinski_Tang', 'Sum_of_Different_Powers',...
    'Sum_Square', 'Zakharov'};
nproblems = size(BoxProblems, 2);
for p = 1:nproblems
    Problem.f = BoxProblems{p};
    [Fmin, Xmin, history] = dDirect(Problem, opts);
    
% Save function evaluations reuslts.    
    DIRECT_filename = sprintf('DIRECT_eval%p', p);
    DIRECT_ID = fopen(DIRECT_filename, 'a');
    fprintf(DIRECT_ID, '%8i\n', history(end, 2));
    fclose(DIRECT_ID);
end
% -------------------------------------------------------------------------
% Test problems with fixed number of variables.
% -------------------------------------------------------------------------
BoxProblems ={'Beale', 'Bohachecsky1', 'Bohachecsky2', 'Bohachecsky3',...
    'Booth', 'Branin', 'Bukin6', 'Colville','Cross_in_Tray',...
    'Drop_wave', 'Easom', 'Eggholder', 'Goldstein_and_Price',...
    'Hartman3', 'Hartman6', 'Holder', 'Hump', 'Langermann', 'Matyas',...
    'McCormick', 'Powell', 'Power_Sum', 'Shekel5', 'Shekel7', 'Shekel10'...
    'Shubert', 'Trid6', 'Trid10'};
nproblems = size(BoxProblems, 2);
for p = 1:nproblems
    Problem.f = BoxProblems{p};
    [Fmin, Xmin, history] = dDirect(Problem, opts);
    
% Save function evaluations reuslts.    
    DIRECT_filename = sprintf('DIRECT_eval%p', p);
    DIRECT_ID = fopen(DIRECT_filename, 'a');
    fprintf(DIRECT_ID, '%8i\n', history(end, 2));
    fclose(DIRECT_ID);
end
% -------------------------------------------------------------------------