% Load options:
% -------------------------------------------------------------------------
opts.maxevals = 10000;   % Maximum number of function evaluations.
opts.maxits   = 10000;   % Maximum number of iterations.
opts.maxdeep  = 10000;   % Maximum number of side divisions.
opts.testflag = 1;       % Testing algorithms with known globalmin vcalue.
opts.showits  = 1;       % Print iteration stats.
opts.tol      = 0.01;    % Allowable relative error if globalmin is set.
% -------------------------------------------------------------------------

% Test problems where the number of variables can be any integer n:
opts.dimension = 2;      % Specify number of variables
% -------------------------------------------------------------------------
BoxProblems ={'Tproblem'};
nproblems = size(BoxProblems, 2);
for p = 1:nproblems
    Problem.f = BoxProblems{p};
    [Fmin, Xmin, history] = dDirect_GLc(Problem, opts);
    
% Save function evaluations reuslts.    
    DIRECT_filename = sprintf('DIRECT_eval%p', p);
    DIRECT_ID = fopen(DIRECT_filename, 'a');
    fprintf(DIRECT_ID, '%8i\n', history(end, 2));
    fclose(DIRECT_ID);
end
% -------------------------------------------------------------------------

% Test problems with fixed number of variables.
% -------------------------------------------------------------------------
BoxProblems ={'Bunnag1', 'Bunnag2', 'Bunnag3', 'Bunnag4', 'Bunnag5',...
    'Bunnag6', 'Bunnag7', 'circle', 'G01', 'G02', 'G04', 'G06', 'G07',...
    'G08', 'G09', 'G10', 'G12', 'G16', 'G18', 'G19', 'G24', 'Genocop9',...
    'Genocop10', 'Genocop11', 'Goldstein_and_Pricec', 'Gomez',... 
    'Himmelblaus', 'Horst1', 'Horst2', 'Horst3', 'Horst4', 'Horst5'...
    ,'Horst6', 'Horst7', 'hs021', 'hs021mod', 'hs024', 'hs035', 'hs036',...
    'hs037', 'hs038', 'hs044', 'hs076', 'P1', 'P2a', 'P2b', 'P2c',...
    'P2d', 'P3a', 'P3b', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10',...
    'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 's224', 's231', 's232',...
    's250', 's251', 's365mod', 'zecevic2', 'zecevic3', 'zecevic4', 'zy2'};
nproblems = size(BoxProblems, 2);
for p = 1:nproblems
    Problem.f = BoxProblems{p};
    [Fmin, Xmin, history] = dDirect_GLc(Problem, opts);
    
% Save function evaluations reuslts.    
    DIRECT_filename = sprintf('DIRECT_eval%p', p);
    DIRECT_ID = fopen(DIRECT_filename, 'a');
    fprintf(DIRECT_ID, '%8i\n', history(end, 2));
    fclose(DIRECT_ID);
end
% -------------------------------------------------------------------------