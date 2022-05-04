# DIRECTGO
**DIRECTGO**: A new **DIRECT**-type `MATLAB` toolbox for derivative-free **G**lobal **O**ptimization [[29]](https://arxiv.org/abs/2107.02205)

---

## Quick overview

The sequential [[1](http://www4.ncsu.edu/~ctk/Finkel_Direct/), [26](https://doi.org/10.1023/A:1019992822938), [27](https://doi.org/10.1016/j.amc.2020.125596)] and parallel implementations [[27]](https://doi.org/10.1016/j.amc.2020.125596) of various DIRECT-type algorithms. The toolbox consists of two main parts:

- **DIRECTGO.mltbx** - `MATLAB` toolbox package containing implementations of DIRECT-type algorithms, including an extensive [DIRECTGOLib](https://github.com/blockchain-group/DIRECTGOLib) [[28]](https://github.com/blockchain-group/DIRECTGOLib) library (maintained separetely) of the box, generally constrained, and practical engineering global optimization problems, often used for benchmarking DIRECT-type algorithms.
- **DIRECTGO.mlappinstall** - A single `MATLAB` app installer ﬁle containing everything necessary to install and run the **DIRECTGO** toolbox, including a graphical user interface (GUI).

Additionally, we provide source files of all implemented algorithms in [Algorithms/](Algorithms/) folder.

## Versions history

- **DIRECTGO** [v1.0.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.0.0) version presented in [[29]](https://arxiv.org/abs/2107.02205) implements 36 different DIRECT-type algorithms.
- Twelve new algorithms, presented in [[30]](https://arxiv.org/abs/2109.14912), are included in **DIRECTGO** [v1.1.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.1.0) (48 algorithms in total)

## Algorithms within DIRECTGO

Classification of 48 implemented DIRECT-type algorithms included in different versions of **DIRECTGO**:


| Version                                                      | Problem type          | Algorithm name & [References]                                |
| ------------------------------------------------------------ | --------------------- | ------------------------------------------------------------ |
| [v1.0.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.0.0) | Box constrained       | DIRECT v4.0 [[1](http://www4.ncsu.edu/~ctk/Finkel_Direct/), [2](https://doi.org/10.1007/BF00941892)], DIRECT-restart [[3]](https://repository.lib.ncsu.edu/handle/1840.4/461), DIRECT-m [[4]]( https://doi.org/10.1007/s10898-006-9029-9), DIRECT-l [[5]](https://doi.org/10.1023/A:1017930332101), DIRECT-rev [[6]](http://www4.ncsu.edu/~ctk/Finkel_Direct/), DIRECT-a [[7]](https://doi.org/10.1007/s10898-012-9952-x), DIRMIN [[8]](https://doi.org/10.1007/s10589-008-9217-2), PLOR [[9]](https://doi.org/10.1007/s10898-015-0364-6), glbSolve [[2](https://doi.org/10.1007/BF00941892), [10](https://www.mat.univie.ac.at/~neum/glopt/mss/BjoeH99.pdf)], glbSolve-sym [[11]](https://doi.org/10.1007/s10898-012-0020-3), glbSolve-sym2 [[11]](https://doi.org/10.1007/s10898-012-0020-3), MrDIRECT [[12]](https://doi.org/10.1007/s10898-016-0447-z), MrDIRECT075 [[13]](https://doi.org/10.1007/s10898-014-0241-8), BIRECT [[14]](https://doi.org/10.1007/s10898-016-0485-6), GB-DISIMPL-C [[15]](https://doi.org/10.1007/s10898-014-0180-4), GB-DISIMPL-V [[15]](https://doi.org/10.1007/s10898-014-0180-4), Gb-BIRECT [[16]](https://doi.org/10.1016/j.eswa.2019.113052), BIRMIN [[16]](https://doi.org/10.1016/j.eswa.2019.113052), Gb-glbSolve [[16]](https://doi.org/10.1016/j.eswa.2019.113052), DISIMPL-C [[17]](https://doi.org/10.1007/s10898-013-0089-3), DISIMPL-V [[17]](https://doi.org/10.1007/s10898-013-0089-3), ADC [[18]](https://doi.org/10.1137/040621132), Aggressive DIRECT [[19]](%5Bdownload%20(psu.edu)%5D(https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.63.280&rep=rep1&type=pdf)), DIRECT-G [[20]](https://doi.org/10.1007/s11590-017-1228-4), DIRECT-L [[20]](https://doi.org/10.1007/s11590-017-1228-4), DIRECT-GL [[20]](https://doi.org/10.1007/s11590-017-1228-4). |
| [v1.1.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.1.0) | Box constrained       | I_DBDP_GL [[30]](https://arxiv.org/abs/2109.14912), I_DBDP_IA [[30]](https://arxiv.org/abs/2109.14912), I_DBDP_IO [[30]](https://arxiv.org/abs/2109.14912), I_DTC_GL [[30]](https://arxiv.org/abs/2109.14912), I_DTC_IA [[30]](https://arxiv.org/abs/2109.14912), I_DTC_IO [[30]](https://arxiv.org/abs/2109.14912), I_DTDV_IA [[30]](https://arxiv.org/abs/2109.14912), I_DTDV_GL [[30]](https://arxiv.org/abs/2109.14912), I_DTDV_IO [[30]](https://arxiv.org/abs/2109.14912), N_DTC_GL [[30]](https://arxiv.org/abs/2109.14912), N_DTC_IA [[30]](https://arxiv.org/abs/2109.14912), N_DTC_IO [[30]](https://arxiv.org/abs/2109.14912). |
| [v1.0.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.0.0) | Linearly constrained  | Lc-DISIMPL-C [[21]](https://doi.org/10.1007/s11590-014-0772-4), Lc-DISIMPL-V [[21]](https://doi.org/10.1007/s11590-014-0772-4). |
| [v1.0.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.0.0) | Generally constrained | DIRECT-L1 [[1]](http://www4.ncsu.edu/~ctk/Finkel_Direct/), DIRECT-GLc [[22]](https://doi.org/10.1007/s00158-018-2181-2), DIRECT-GLce [[22]](https://doi.org/10.1007/s00158-018-2181-2), DIRECT-GLce-min [[22]](https://doi.org/10.1007/s00158-018-2181-2). |
| [v1.0.0](https://github.com/blockchain-group/DIRECTGO/releases/tag/v1.0.0) | Hidden constraints    | DIRECT-NAS [[1]](http://www4.ncsu.edu/~ctk/Finkel_Direct/), DIRECT-Barrier [[23]](%5BModification%20of%20the%20DIRECT%20Algorithm%20(ncsu.edu)%5D(https://repository.lib.ncsu.edu/bitstream/handle/1840.16/3920/etd.pdf?sequence=1)), subDIRECT-Barrier [[24]](https://doi.org/10.1016/j.energy.2017.03.047), DIRECT-GLh [[25]](https://doi.org/10.1007/s11590-021-01726-z). |

## Quick user guide

After installation of the `MATLAB` toolbox (using **DIRECTGO.mltbx**), all implemented DIRECT-type algorithms and test problems can be freely accessed in the command window of `MATLAB`. Unlike using GUI, algorithms from the command line require more programming knowledge, and configurations must be done manually. All algorithms can be run using the same style and syntax:

```matlab
1. f_min = algorithm(P);
2. f_min = algorithm(P, OPTS);
3. f_min = algorithm(P, OPTS, D);
4. [f_min, x_min] = algorithm(P, OPTS, D);
5. [f_min, x_min, history] = algorithm(P, OPTS, D);
```

The left side of the equations specifies the output parameters. After the termination, the algorithm returns the best objective value `f_min`, solution point `x_min`, and history of the algorithmic performance during all iterations `history`. The information presented here is: the iteration number, the total number of objective function evaluations, the current minimal value, and the total execution time.

On the right side, the algorithm name `algorithm` and at least one input parameters are needed to specify. The first one is the problem structure `P` consisting of an objective function:

```matlab
>> P.f = 'objfun';
```

If the problem involves additional constraints, they also must be specified:

```matlab
>> P.constraint = 'confun';
```

The second parameter is an optional variable `OPTS` to customize the default settings. The last parameter is the bound constraints for each dimension/variable:

```matlab
D (i,1) ≤ x_i ≤ D (i,2), i = 1...n;
```

The latter `D` parameter is necessary for the algorithm, but can be passed to it in other ways.


### Example of box constrained global optimization algorithm usage

---

Any DIRECT-type algorithm in the **DIRECTGO** toolbox can be called in `MATLAB` using the same sequences, which are discussed earlier. Let’s consider the example using the *PLOR* algorithm, solving **Bukin6** test function. Figure 1. (left side) shows the **Bukin6** test function plot over its domain. The function should be defined as:

```matlab
function y = Bukin6(x)                 % Extract info from the function
    if nargin == 0
        y.nx = 2;                      % Dimension of the problem
        xl = [-15; -3];
        y.xl = @(i) xl(i);             % Lower bounds for each variable
        xu = [5; 3];
        y.xu = @(i) xu(i);     	       % Upper bounds for each variable
	      y.fmin = @(nx) get_fmin(nx);   % Known solution value
        y.xmin = @(nx) get_xmin(nx);   % Known solution point
        return
    end
    if size(x, 2) > size(x, 1)         % If x is a row transpose to column 
        x = x'; 
    end
    term1 = 100*sqrt(abs(x(2) - 0.01*x(1)^2));
    term2 = 0.01*abs(x(1) + 10);
    y = term1 + term2;                 % Return function value at x
end

function fmin = get_fmin(~)
    fmin = 0;
end

function xmin = get_xmin(~)
    xmin = [-10; 1];
end
```

Each test problem in the **DIRECTGOLib** stores the information about the problem structure together with the objective function. In this case, in the `Bukin6.m` file, the following information is stored: i) the dimensionality of the problem; ii) the lower and upper bounds for each variable; iii) the objective function value of the known solution; iv) the solution point.
The optimization problem is passed to the algorithm as part of a `P` structure. For a simple, box-constrained problem like this, only one field of the `P` structure is needed:
For some problems, the optimum might depend on the number of variables, therefore the solution values and points are returned as functions for all test problems in **DIRECTGOLib**.

```matlab
» P.f = 'Bukin6';
```

When a user wants to perform calculations using different from the default settings, the `OPTS` structure should be used:

```matlab
>> opts.maxevals = 50;       % Maximal number of function evaluations
>> opts.maxits = 100;        % Maximal number of iterations
>> opts.testflag = 1;        % 1 if global minima is known, 0 otherwise
>> opts.tol = 0.01;          % Tolerance for termination if testflag = 1
```

Now we are ready to call the dynamic data structure based *PLOR* implementation `dPlor.m` to solve this problem:

```matlab
>> [f_min, x_min, history] = dPlor(P, OPTS);
```

The iterative output stored in the history parameter contains the following information:

```matlab
>> history
        history =
            1.0000        5.0000       16.7833       0.0023
            2.0000        7.0000       16.7833       0.0030
            3.0000       13.0000        5.6500       0.0039
            4.0000       19.0000        5.6500       0.0046
            5.0000       27.0000        1.9537       0.0053
            6.0000       33.0000        1.9537       0.0060
            7.0000       41.0000        0.7167       0.0070
            8.0000       47.0000        0.7167       0.0077
            9.0000       55.0000        0.3060       0.0086
```

Here, the first column shows the iteration number, while the second is the total number of function evaluations. The third column shows how the best objective function value improves at each iteration, while the last column shows the execution time in seconds. The *PLOR* algorithm was terminated when the maximum number of function evaluations (opts.maxevals = 50) exceeded.
The convergence plot is shown on the right side of Fig. 1, while the left side illustrates the **Bukin6** test function over its domain.

![bukin6](img/bukin6.png)

### Example of constrained global optimization algorithm usage

---

Any DIRECT-type algorithmic implementation for constrained global optimization problems can be used, following the same principle presented earlier. For constrained problems, implemented algorithms extract additional information from the functions, such as: the number of inequality constraints; the number of equality constraints; and the constraint functions. Let us take the **G06** problem as an example, which is defined in the following way:

```matlab
function y = G06(x)                           
    if nargin == 0			                  % Extract info from the function
        y.nx = 2;                             % Dimension of the problem
        y.ng = 2;                             % Number of g(x) constraints
        y.nh = 0;                             % Number of h(x) constraints
        xl = [13, 0];
        y.xl = @(i) xl(i);                    % Lower bounds for each variable
        y.xu = @(i) 100;		                  % Upper bounds for each variable
        y.fmin = @(nx) get_fmin(nx);          % Known solution value
        y.xmin = @(nx) get_xmin(nx);          % Known solution point
        y.confun = @(i) G06c(i);              % Constraint functions
        return
    end
    if size(x, 2) > size(x, 1)                % If x is a row transpose to column 
        x = x'; 
    end
    y = (x(1) - 10)^3 + (x(2) - 20)^3;        % Return function value at x
end

function [c, ceq] = G06c( x )
    c(1) = -(x(1) - 5)^2 - (x(2) - 5)^2 + 100;
    c(2) = (x(1) - 6)^2 + (x(2) - 5)^2 - 82.81;
    ceq = [];
end

function fmin = get_fmin(~)
    fmin = -6961.8138751273809248;
end

function xmin = get_xmin(~)
    xmin = [14.0950000002011322; 0.8429607896175201];
end
```

Same as in the previous example the optimization problem is passed to the algorithm as part of a `P` structure. Once again, only one field of the `P` structure is needed:

```matlab
>> P.f = 'G06';
```

Next, assume that the user wants to stop the search as soon as the known solution is found to be within 0.01% error. The structure of `OPTS` in the `MATLAB` command window should be given as follows:

```matlab
>> opts.maxevals = 10000;  % Maximal number of function evaluations
>> opts.maxits = 1000;     % Maximal number of iterations
>> opts.testflag = 1;      % 1 if global minima is known, 0 otherwise
>> opts.tol = 0.01;        % Tolerance for termination if testflag = 1
>> opts.showits = 1;       % Print iteration status
```

The desired algorithm can then be run using the following sequence:

```matlab
>> [f_min, x_min, history] = dDirect_GLc(P, OPTS);
```

Since `opts.showits` has been set to 1 in the `OPTS` structure, the status of each iteration will be printed in the `MATLAB` command window. The *DIRECT-GLc* algorithm used in this example works in two phases. Since the initial sampling points of the `G06` problem do not satisfy the constraints using *DIRECT-GLc*, the algorithm switches to the second phase, which is indicated by the first printed line. In the latter phase, the algorithm searches for at least one point where the constraints are satisfied. When the algorithm finds such a point, it switches to phase one and continues trying to find a better solution using the auxiliary function based approach. The last line prints the reason for stopping the algorithm. Since the solution was found with a 0.01% error, the algorithm stopped at iteration 14 after finding the solution 𝑓min with the value −6901.5099081387.

```
Phase_II - searching feasible point:
con viol: 2404.4400000000 fn evals: 5
con viol: 515.5511111111  fn evals: 7
...
con viol: 0.1374240038    fn evals: 123
con viol: 0.0000000000    fn evals: 159  f_min: -5612.1483164940
Phase_I - Improve feasible solution:
Iter: 1     f_min: -5886.5625227848     time(s): 0.05935     fn evals: 197
Iter: 2     f_min: -5931.8554991123     time(s): 0.06473     fn evals: 241
...
Iter: 13    f_min: -6873.0583159376     time(s): 0.13197     fn evals: 947
Iter: 14    f_min: -6901.5099081387     time(s): 0.13869     fn evals: 1027
Minima was found with Tolerance: 1
```

Let us consider a different way in which a user can use DIRECT-type algorithms. If the user wants to solve new problems that are not available in **DIRECTGOLib** without defining the functions as described in earlier, this can be done in other ways. Let us take the same `G06` problem as an example. To the structure `P.f`, the name of an m-file which should compute the value of the `G06` objective function must be given.

```matlab
function y = G06(x)
    y = (x(1) - 10)^3 + (x(2) - 20)^3;
end
```

To the structure `P.constraint`, the name of an m-file which should compute the vector of the **G06** constrain functions must be given.

```matlab
function [c, ceq] = G06c(x)
    c(1) = -(x(1) - 5)^2 - (x(2) - 5)^2 + 100;
    c(2) = (x(1) - 6)^2 + (x(2) - 5)^2 - 82.81;
    ceq = [];
end
```

After the objective and constraint functions are passed to the structure `P`:

```matlab
>> P.f = 'G06';
>> P.constraint = 'G06c';
```

The next necessary parameter to be passed is the optimization domain `D`. `G06` is a second dimension test problem, and domain `D` should look like a 2x2 matrix:

```matlab
>> D = [13, 100; 0, 100];
```

The first column should indicate the lower bounds for `x`, and the second column the upper bounds. Next, assume that we going to use the same `OPTS` structure, which we already defined earlier. However, for the algorithm to find a solution with the desired 0.01% error, it needs to specify the solution of the `G06` problem in the `OPTS` structure.

```matlab
>> opts.globalmin = -6961.81387512; % Known global solution value
```

Then, the desired algorithm can be used with the following sequence:

```matlab
>> [f_min, x_min, history] = dDirect_GLc(P, OPTS, D);
```

This results in the same iterative sequence of the algorithm as in the previous example in this section. The *DIRECT-GLc* algorithm will terminate after 14 iterations and 1027 evaluations of the objective function.

**Parallel algorithm usage**

This section briefly explains how to use parallel versions of the algorithms. Assume a user wishes to use parallel code for the *PLOR* algorithm. First, a parallel implementation of the *PLOR* algorithm `parallel_dPlor.m` should be chosen. Next, a user should specify the number of workers (computational threads). For parallel *PLOR*, it is reasonable to select 2, as only two potential optimal hyper-rectangles are selected per iteration. In this case, `MATLAB` parallel pool size should be specified using the `parpool` command, after which the parallel algorithm should be executed:

```matlab
>> parpool(2);
>> [f_min, x_min, history] = parallel_dPlor(P, OPTS);
```

By default, the `parpool` command starts the `MATLAB` pool on the local machine with one worker per physical CPU core. Using `parpool(2)`, we limit the number of workers to 2. After this, the parallel code is executed using both workers. However, it should be taken into account that creating parallel `parpool` takes some time. Therefore, using the parallel *PLOR* algorithm is inefficient in solving simple problems. The use of parallel codes should address higher-dimensionality, more expensive optimization problems. When all necessary calculations in parallel mode are finished, the following command:

```matlab
>> delete(gcp);
```

shuts down the parallel pool.

## Using Scripts

### Reproducing results from [[29]](https://doi.org/10.48550/arXiv.2107.02205)

Four scripts in the folder [Scripts/TOMS](https://github.com/blockchain-group/DIRECTGO/tree/main/Scripts/TOMS) can be used to reproduce results presented in the manuscript: [DIRECTGO: A new DIRECT-type MATLAB toolbox for derivative-free global optimization](https://doi.org/10.48550/arXiv.2107.02205), and given in [Results/TOMS](https://github.com/blockchain-group/DIRECTGO/tree/main/Results/TOMS) folder. The scripts automatically download the required version of [DIRECTGOLib](https://github.com/blockchain-group/DIRECTGOLib) for experiments.

- `SolveBoxProblems.m` - can be used to repeat experiments for box-constrained problems presented in TABLE 3 and 4 [[29]](https://doi.org/10.48550/arXiv.2107.02205);
- `SolveGeneralProblems.m` - can be used to repeat experiments for box-constrained problems presented in TABLE 5 [[29]](https://doi.org/10.48550/arXiv.2107.02205);
- `SolveGeneralPracticalProblems.m` - can be used to repeat experiments for box-constrained problems presented in TABLE 6-10 [[29]](https://doi.org/10.48550/arXiv.2107.02205);
- `SolveBoxPracticalProblems.m` - can be used to repeat experiments for box-constrained problems presented in TABLE 11 and 12 [[29]](https://doi.org/10.48550/arXiv.2107.02205).

### Reproducing results from [[30]](https://arxiv.org/abs/2109.14912)

The script in the folder [Scripts/JOGO](https://github.com/blockchain-group/DIRECTGO/tree/main/Scripts/JOGO) can be used to reproduce results presented in the manuscript: [An empirical study of various candidate selection and partitioning techniques in the DIRECT framework](https://arxiv.org/abs/2109.14912) and given in [Results/JOGO](https://github.com/blockchain-group/DIRECTGO/tree/main/Results/JOGO) folder. The scripts automatically download the required version of [DIRECTGOLib](https://github.com/blockchain-group/DIRECTGOLib) for experiments.

- `SolveDIRECTGOlib.m` - can be used to repeat experiments presented in TABLE 2 [[30]](https://arxiv.org/abs/2109.14912).

### Reproducing HALRECT results

The script in the folder [Scripts/COA](https://github.com/blockchain-group/DIRECTGO/tree/main/Scripts/COA) can be used to reproduce results presented in the manuscript: *"Lipschitz-inspired HALRECT Algorithm for Derivative-free Global Optimization"* and given in [Results/COA](https://github.com/blockchain-group/DIRECTGO/tree/main/Results/COA) folder. The scripts automatically download the required version of [DIRECTGOLib](https://github.com/blockchain-group/DIRECTGOLib) for experiments.

- `SolveHALRECT.m` - can be used to repeat experiments presented in TABLE 3, 4 and FIGURE 10.

## Citing DIRECTGO

Please use the following bibtex entries, if you consider to cite `DIRECTGO` toolbox:

```latex
@misc{Stripinis2022:DirectGOv1.0,
  title        = {{DIRECTGO: A new DIRECT-type MATLAB toolbox for derivative-free global optimization}}, 
  author       = {Linas Stripinis and Remigijus Paulavi{\v c}ius},
  year         = {2022},
  publisher    = {GitHub},
  version      = {v1.0},
  howpublished = {\url{https://github.com/blockchain-group/DIRECTGO}},
}

@misc{Stripinis2022:dgo,
	title     = {{DIRECTGO: A new DIRECT-type MATLAB toolbox for derivative-free global optimization}}, 
	author    = {Linas Stripinis and Remigijus Paulavi{\v{c}}ius},
	year      = {2022},
	eprint    = {arXiv:2107.02205v2},
	publisher = {arXiv},
	url       = {https://arxiv.org/abs/2107.02205}
}
```

## References

1. D. E. Finkel. 2004. MATLAB source code for DIRECT. http://www4.ncsu.edu/~ctk/Finkel_Direct/. Online; accessed: 2017-03-22.
2. D. R. Jones, C. D. Perttunen, and B. E. Stuckman. 1993. Lipschitzian Optimization Without the Lipschitz Constant. Journal of Optimization Theory and Application 79, 1 (1993), 157–181. https://doi.org/10.1007/BF00941892
3. D. Finkel and C. Kelley. 2004. An Adaptive Restart Implementation of DIRECT. [An adaptive restart implementation of direct (ncsu.edu)](https://repository.lib.ncsu.edu/handle/1840.4/461) In Technical report CRSC-TR04-30. Center for Research in Scientific Computation, North Carolina State University, Raleigh, 1–16.
4. D. E. Finkel and C. T. Kelley. 2006. Additive scaling and the DIRECT algorithm. Journal of Global Optimization 36, 4 (2006), 597–608. [https://doi.org/10.1007/s10898-006-9029-9]( https://doi.org/10.1007/s10898-006-9029-9)
5. J. M. Gablonsky and C. T. Kelley. 2001. A locally-biased form of the DIRECT algorithm. Journal of Global Optimization 21, 1 (2001), 27–37. https://doi.org/10.1023/A:1017930332101
6. D. R. Jones. 2001. The Direct Global Optimization Algorithm. http://www4.ncsu.edu/~ctk/Finkel_Direct/. In The Encyclopedia of Optimization, Christodoulos A. Floudas and Panos M. Pardalos (Eds.). Kluwer Academic Publishers, Dordrect, 431–440.
7. Qunfeng Liu. 2013. Linear scaling and the DIRECT algorithm. Journal of Global Optimization 56 (2013), 1233–1245. Issue 3. https://doi.org/10.1007/s10898-012-9952-x
8. G. Liuzzi, S. Lucidi, and V. Piccialli. 2010. A DIRECT-based approach exploiting local minimizations for the solution of large-scale global optimization problems. Computational Optimization and Applications 45 (2010), 353–375. Issue 2. https://doi.org/10.1007/s10589-008-9217-2
9. Jonas Mockus, Remigijus Paulavičius, Dainius Rusakevičius, Dmitrij Šešok, and Julius Žilinskas. 2017. Application of Reduced-set Pareto-Lipschitzian Optimization to truss optimization. Journal of Global Optimization 67, 1-2 (2017), 425–450. https://doi.org/10.1007/s10898-015-0364-6
10. Mattias Björkman and Kenneth Holmström. 1999. Global Optimization Using the DIRECT Algorithm in Matlab. https://www.mat.univie.ac.at/~neum/glopt/mss/BjoeH99.pdf. Advanced Modeling and Optimization 1, 2 (1999), 17–37.
11. Ratko Grbić, Emmanuel Karlo Nyarko, and Rudolf Scitovski. 2013. A modification of the direct method for Lipschitz global optimizatio n for a symmetric function. Journal of Global Optimization 57, 4 (2013), 1193–1212. https://doi.org/10.1007/s10898-012-0020-3
12. Qunfeng Liu, Guang Yang, Zhongzhi Zhang, and Jinping Zeng. 2017. Improving the convergence rate of the DIRECT global optimization algorithm. Journal of Global Optimization 67, 4 (2017), 851–872. https://doi.org/10.1007/s10898-016-0447-z 
13. Qunfeng Liu, Jinping Zeng, and Gang Yang. 2015. MrDIRECT: a multilevel robust DIRECT algorithm for global optimization problems. Journal of Global Optimization 62, 2 (2015), 205–227. https://doi.org/10.1007/s10898-014-0241-8
14. Remigijus Paulavičius, Lakhdar Chiter, and Julius Žilinskas. 2018. Global optimization based on bisection of rectangles, function values at diagonals, and a set of Lipschitz constants. Journal of Global Optimization 71, 1 (2018), 5–20. https://doi.org/10.1007/s10898-016-0485-6
15. Remigijus Paulavičius, Ya. D. Sergeyev, Dmitri E. Kvasov, and Julius Žilinskas. 2014. Globally-biased DISIMPL algorithm for expensive global optimization. Journal of Global Optimization 59, 2-3 (2014), 545–567. https://doi.org/10.1007/s10898-014-0180-4
16. Remigijus Paulavičius, Ya. D. Sergeyev, Dmitri E. Kvasov, and Julius Žilinskas. 2020. Globally-biased BIRECT algorithm with local accelerators for expensive global optimization. Expert Systems with Applications 144 (2020), 11305. https://doi.org/10.1016/j.eswa.2019.113052
17. Remigijus Paulavičius and Julius Žilinskas. 2013. Simplicial Lipschitz optimization without the Lipschitz constant. Journal of Global Optimization 59, 1 (2013), 23–40. https://doi.org/10.1007/s10898-013-0089-3
18. Ya. D. Sergeyev and Dmitri E. Kvasov. 2006. Global search based on diagonal partitions and a set of Lipschitz constants. SIAM Journal on Optimization 16, 3 (2006), 910–937. https://doi.org/10.1137/040621132
19. Chuck A. Baker, Layne T. Watson, Bernard Grossman, William H. Mason, and Raphael T. Haftka. 2001. [Parallel Global Aircraft Configuration Design Space Exploration]([download (psu.edu)](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.63.280&rep=rep1&type=pdf)). Nova Science Publishers, Inc., USA, 79–96.
20. Linas Stripinis, Remigijus Paulavičius, and Julius Žilinskas. 2018. Improved scheme for selection of potentially optimal hyper-rectangles in DIRECT. Optimization Letters 12, 7 (2018), 1699–1712. https://doi.org/10.1007/s11590-017-1228-4
21. Remigijus Paulavičius and Julius Žilinskas. 2016. Advantages of simplicial partitioning for Lipschitz optimization problems with linear constraints. Optimization Letters 10, 2 (2016), 237–246. https://doi.org/10.1007/s11590-014-0772-4
22. Linas Stripinis, Remigijus Paulavičius, and Julius Žilinskas. 2019. Penalty functions and two-step selection procedure based DIRECT-type algorithm for constrained global optimization. Structural and Multidisciplinary Optimization 59, 6 (2019), 2155–2175. https://doi.org/10.1007/s00158-018-2181-2
23. J. M. Gablonsky. 2001. [Modifications of the DIRECT Algorithm]([Modification of the DIRECT Algorithm (ncsu.edu)](https://repository.lib.ncsu.edu/bitstream/handle/1840.16/3920/etd.pdf?sequence=1)). Ph.D. Dissertation. North Carolina State University.
24. Jonggeol Na, Youngsub Lim, and Chonghun Han. 2017. A modified DIRECT algorithm for hidden constraints in an LNG process optimization. Energy 126 (2017), 488–500. https://doi.org/10.1016/j.energy.2017.03.047
25. Linas Stripinis and Remigijus Paulavičius. 2021. A new DIRECT-GLh algorithm for global optimization with hidden constraints. Optimization Letters 15, 6 (2021), 1865–1884. https://doi.org/10.1007/s11590-021-01726-z
26. Jian He, Layne T. Watson, Naren Ramakrishnan, Clifford A. Shaffer, Alex Verstak, Jing Jiang, Kyung Bae, and William H. Tranter. 2002. Dynamic Data Structures for a DIRECT Search Algorithm. Computational Optimization and Applications 23, 1 (2002), 5–25. https://doi.org/10.1023/A:1019992822938
27. Linas Stripinis, Julius Žilinskas, Leocadio G. Casado, and Remigijus Paulavičius. 2021. On MATLAB experience in accelerating DIRECT-GLce algorithm for constrained global optimization through dynamic data structures and parallelization. Appl. Math. Comput. 390 (2021), 1–17. https://doi.org/10.1016/j.amc.2020.125596
28. Linas Stripinis and Remigijus Paulavičius. 2022. DIRECTGOLib - DIRECT Global Optimization test problems Library. https://doi.org/10.5281/zenodo.6491863
29. Stripinis, L., Paulavičius, R.: DIRECTGO: A new DIRECT-type MATLAB toolbox for derivative-free global optimization (2022). URL https://arxiv.org/abs/2107.02205
30. Stripinis, L., Paulavičius, R.: An empirical study of various candidate selection and partitioning techniques in the DIRECT framework. arXiv (2021). https://doi.org/10.48550/ARXIV.2109.14912. https://arxiv.org/abs/2109.14912
