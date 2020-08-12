function [fmin,x] = DIRECT_Benchmark_P7(options)
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

% 1. Establish bounds for variables
bounds = [0 1;0 1;0 1;0 1];

% 2. Send options to Direct

% 2a. NEW!
% Pass Function as part of a Matlab Structure
cost = [0   0   0   1   5];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 4;

g1 = [0    0    0     0    1.4000;
    1.0000    0    0          0   -5.0000;
         0   0    0      1.0000   -1.2500];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = 1;
  
g2 = [0   0    0      0   -1.4000;
    1.0000    0    0     0    5.0000;
         0  0    0       1.0000   -1.2500];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = 1;

g3 = [0      0     0    0    1.5000;
    0  1.0000    0     0   -5.0000;
    0     0  0  1.0000   -1.0000];
Problem.constraint(3).func = @(x)evaluate_function(g3,x);
Problem.constraint(3).penalty = 1;

g4 = [   0      0     0    0   -1.5000;
    0   1.0000     0    0    5.0000;
      0   0  0  1.0000   -1.0000];
Problem.constraint(4).func = @(x)evaluate_function(g4,x);
Problem.constraint(4).penalty = 1;

g5 = [ 0    0     0         0    0.8000;
    0    0   1.0000         0   -5.0000;
      0    0      0    1.0000   -1.0000];
Problem.constraint(5).func = @(x)evaluate_function(g5,x);
Problem.constraint(5).penalty = 1;

g6 = [   0    0     0         0   -0.8000;
    0    0    1.0000         0    5.0000;
     0    0    0    1.0000   -1.0000];
Problem.constraint(6).func = @(x)evaluate_function(g6,x);
Problem.constraint(6).penalty = 1;

g7 = [4           0           0        0        -625;
    4           4           0        0      390625;
    0           4           1        0       -3125];
Problem.constraint(7).func = @(x)evaluate_function(g7,x);
Problem.constraint(7).penalty = 1;

g8 = [4           0           0        0        625;
    4           4           0        0      -390625;
    0           4           1        0       3125];
Problem.constraint(8).func = @(x)evaluate_function(g8,x);
Problem.constraint(8).penalty = 1;

% 3. Call DIRECT
[fmin,x,~] = Direct(Problem,bounds,options);

end

