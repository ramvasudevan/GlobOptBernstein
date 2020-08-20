function [fmin,x,time_spent,feasibility] = DIRECT_Benchmark_P2(options)
% Benchmark problem setup for DIRECT

%------------------------------------------------------%
% DIRECT Function Call                                 %
%                                                      %
% You will need the accompanying programs              %
% Direct.m, and gp_con.m to run this code              %
% successfully.                                        %
%                                                      %
% These codes can be found at                          %
% http://www4.ncsu.edu/~definkel/research/index.html   %
%------------------------------------------------------%

% 0. create default options if necessary
if nargin < 1
    options.maxevals  = 50000;
    options.maxits    = 46; % TUNED BY HAND
    options.testflag  = 0;
    options.showits   = 0;
    options.tol = 0.1369 ;
end

% 1. Establish bounds for variables
bounds = [0 1;0 1];
constraint_penalty = 1e5 ; % TUNED BY HAND

% 2. Send options to Direct

% 2a. NEW!
% Pass Function as part of a Matlab Structure
% Problem.f = 'gp_con';
cost = [0, 0, -7973;
        1, 0, 2349;
        2, 0, 68121;
        3, 0, 658500;
        0, 1, 120000;
        0, 2, -600000;
        0, 3, 1000000];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 2;

g1 = [0, 0, 11;
      1, 0, -1392;
      2, 0, -7569;
      0, 1, 1000;
      0, 2, -10000];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = constraint_penalty ;
  
g2 = [0, 0, -8.81;
      1, 0, 1218;
      2, 0, 7569;
      0, 1, -1000;
      0, 2, 10000];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = constraint_penalty ;

% 3. Call DIRECT
time_spent_start = tic ;
[fmin,x,~] = Direct(Problem,bounds,options);
time_spent = toc(time_spent_start) ;

% 4. Evaluate constraints on output
feasibility = true ;
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

