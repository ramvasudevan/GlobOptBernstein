function [fmin,x] = DIRECT_Benchmark_P5(options)
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
bounds = [0 1;0 1;0 1];

% 2. Send options to Direct

% 2a. NEW!
% Pass Function as part of a Matlab Structure
cost = [0,  0,  0,  -5;
        0,  0,  1,   10];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 4;

g1 = [     0           0           0        -149;
           1           0           0        2180;
           2           0           0       -5800;
           3           0           0        4000;
           0           1           0        -200;
           1           1           0         400;
           0           0           1         -10];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = 1;
  
g2 = [     0           0           0         159;
           1           0           0       -2180;
           2           0           0        5800;
           3           0           0       -4000;
           0           1           0         200;
           1           1           0        -400;
           0           0           1         -10];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = 1;

g3 = [     0           0           0        -237
           1           0           0        -400
           2           0           0         200
           0           1           0        2540
           0           2           0       -6000
           0           3           0        4000
           1           1           0         400
           0           0           1         -10];
Problem.constraint(3).func = @(x)evaluate_function(g3,x);
Problem.constraint(3).penalty = 1;

g4 = [     0           0           0         247
           1           0           0         400
           2           0           0        -200
           0           1           0       -2540
           0           2           0        6000
           0           3           0       -4000
           1           1           0        -400
           0           0           1         -10];
Problem.constraint(4).func = @(x)evaluate_function(g4,x);
Problem.constraint(4).penalty = 1;

% 3. Call DIRECT
[fmin,x,~] = Direct(Problem,bounds,options);

end

