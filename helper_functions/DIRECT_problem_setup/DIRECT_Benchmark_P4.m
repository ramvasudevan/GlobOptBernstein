function [fmin,x] = DIRECT_Benchmark_P4(options)
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
    options.maxits    = 40;
    options.testflag  = 0;
    options.showits   = 1;
    options.tol = 1.7e-06 ;
end

% 1. Establish bounds for variables
bounds = [0 1;0 1;0 1];
constraint_penalty = 1e3 ;

% 2. Send options to Direct

% 2a. NEW!
% Pass Function as part of a Matlab Structure
cost = [1     0     0    -4;
        0     1     0    10;
        0     0     1    -3];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 3;

g1 = [0     0     0    -4;
      1     0     0     2;
      0     1     0    10;
      0     0     1     3];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = constraint_penalty;
  
g2 = [0     0     0    -6;
      0     1     0    30;
      0     0     1     3];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = constraint_penalty;

g3 = [0     0     0   -24;
     1     0     0    40;
     2     0     0   -16;
     0     1     0   -90;
     0     2     0  -200;
     1     1     0    80;
     0     0     1    39;
     0     0     2   -18;
     1     0     1   -24;
     0     1     1    60];
Problem.constraint(3).func = @(x)evaluate_function(g3,x);
Problem.constraint(3).penalty = constraint_penalty;

% 3. Call DIRECT
[fmin,x,~] = Direct(Problem,bounds,options);

end

