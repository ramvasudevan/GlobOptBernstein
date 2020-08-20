function [fmin,x,time_spent,feasibility] = DIRECT_benchmark_helper_function(problem_index,show_iterations,problem_accuracy)
% [fmin,x,time_spent,feasibility = DIRECT_benchmark_helper_function(problem_index)
%
% Given the cost, constraints, and problem index, call DIRECT with the
% appropriate options. This function returns the optimal value found, the
% optimizer x, the amount of time spent, and whether or not DIRECT found a
% feasible answer
%
% NOTE this uses a slightly modified DIRECT function (which returns the
% best answer found, as opposed to failing inelegantly)
%
% Author: Shreyas Kousik & Bohao Zhang
% Created: 20 Aug 2020
% Updated: nah

% set up solver options
options.maxevals  = 50000;
options.testflag  = 0;

if nargin > 1
    options.showits = show_iterations ; % set to 1 for tuning direct's hyperparams
else
    options.showits = 0 ;
end

switch problem_index
    % NOTE these maxits and constraint penalty values have been tuned by
    % hand, since DIRECT has to make a horrible tradeoff between cost and
    % constraints :(
    case 1
        options.maxits = 63 ;
        options.tol = 7.0e-7 ;
        constraint_penalty = 2.075 ;
    case 2
        options.maxits = 46;
        options.tol = 0.1369 ;
        constraint_penalty = 1e5 ;
    case 3
        options.maxits = 35 ;
        options.tol = 2e-06 ;
        constraint_penalty = 1e8 ;
    case 4
        options.maxits = 45 ;
        options.tol = 1.7e-06 ;
        constraint_penalty = 3 ;
    case 5
        options.maxits = 22 ;
        options.tol = 1.7e-06 ;
        constraint_penalty = 0.51 ;
    case 6
        options.maxits = 35 ;
        options.tol = 4.2677e-04 ;
        constraint_penalty = 1 ; % this works surprisingly well!
    case 7
        options.maxits = 40 ;
        options.tol = 5e-7 ;
        constraint_penalty = 1e1 ;
    case 8
        options.maxits = 40 ;
        options.tol = 1.3445e-4 ;
        constraint_penalty = 1 ;
    otherwise
        error('Please choose a problem index between 1 and 8')
end

% set accuracy if given by user
if nargin > 2
    disp('Using provided tolerance') 
    options.tol = problem_accuracy ;
else
    disp('Using default tolerance')
end

% create DIRECT problem structure
eval(strcat('[cost,constraints] = setup_problem_matrix_P',num2str(problem_index),'();')) ;
Problem.f = @(x) evaluate_function(cost,x) ;
Problem.numconstraints = length(constraints) ;
for idx = 1:Problem.numconstraints
    Problem.constraint(idx).func = @(x) evaluate_function(constraints{idx},x) ;
    Problem.constraint(idx).penalty = constraint_penalty ;
end

problem_dimension = size(cost,2) - 1 ;
bounds = repmat([0 1],problem_dimension,1) ;

% call DIRECT
time_spent_start = tic ;
[fmin,x,~] = Direct(Problem,bounds,options);
time_spent = toc(time_spent_start) ;

% evaluate constraints on output
feasibility = true ; % optimism!
for idx = 1:Problem.numconstraints
    x_feas_idx = Problem.constraint(idx).func(x) <= 0 ;
    feasibility = feasibility && x_feas_idx ;
end
if feasibility
    disp('DIRECT converged to a feasible solution!')
else
    disp('DIRECT converged to an infeasible solution!')
end
end