% Load path:
% -------------------------------------------------------------------------
parts = strsplit(pwd, filesep);
parts{end} = 'Algorithms';
parent_path = strjoin(parts(1:end), filesep);
addpath(parent_path); 
parts{end} = 'DIRECTlib';
parts{end + 1} = 'Practical';
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
% Practical problems with box constraints.
% -------------------------------------------------------------------------
BoxProblems ={'Non_Reg_n3_T10', 'Non_Reg_n3_T100', 'Non_Reg_n6_T10',...
    'Non_Reg_n6_T100', 'Non_Reg_n9_T10', 'Non_Reg_n9_T100'};
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
% Practical problems with general constraints.
% -------------------------------------------------------------------------
BoxProblems ={'NASA_speed_reducer', 'Pressure_Vessel',...
    'Tension_Compression_Spring', 'Three_bar_truss', 'Welded_Beam'};
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