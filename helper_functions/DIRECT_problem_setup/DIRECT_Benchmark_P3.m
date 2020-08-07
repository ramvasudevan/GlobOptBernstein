function [fmin,x] = DIRECT_Benchmark_P3(options)
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
bounds = [0 1;0 1];

% 2. Send options to Direct

% 2a. NEW!
% Pass Function as part of a Matlab Structure
% Problem.f = 'gp_con';
cost = [0,   0,   -10;
        1,   0,    20];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 2;

g1 = [0,     0,   110;
     1,     0,  -400;
     2,     0,   400;
     0,     1,   -20];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = 1;
  
g2 = [[   0,   0,                              119]
        [ 1.0,   0,                           -680]
        [ 2.0,   0,                           1280]
        [ 3.0,   0,                           -800]
        [   0, 1.0,                              2]];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = 1;

% 3. Call DIRECT
[fmin,x,~] = Direct(Problem,bounds,options);

end

