%------------------------------------------------------%
% DIRECT Sample Call                                   %
% Purpose: Find the Global Minimum value of            %
%          the Goldstein-Price Test Function           %
%          using the DIRECT algorithm                  %
%          subject to the constraint                   %
%          x^2 + y^2 - 1 <= 0                          %
%          USES NAS STRATEGY                           %
%                                                      %
% You will need the accompanying programs              %
% Direct.m, and gp_con.m to run this code              %
% successfully.                                        %
%                                                      %
% These codes can be found at                          %
% http://www4.ncsu.edu/~definkel/research/index.html   %
%------------------------------------------------------%
clear all;
% 1. Establish bounds for variables
% bounds = [-2 2;-2 2];
bounds = [0 1;0 1];

% 2. Send options to Direct
options.maxevals  = 100;
options.testflag  = 0;
options.globalmin = -5.5080132636;
options.showits   = 1;
options.tol       = 0.01;

% NEW Update options to use implicit constraints
options.impcons = 0;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
% Problem.f = 'gp_con';
problem_data = [1,     0,    -3;
                0,     1,    -4];
Problem.f = @(x)evaluate_function(problem_data,x);
g1 = [0,     0,    -2;
     2,     0,   -72;
     3,     0,   216;
     4,     0,  -162;
     0,     1,     4];
  
g2 = [0,     0,   -36;
     1,     0,   288;
     2,     0,  -792;
     3,     0,   864;
     4,     0,  -324;
     0,     1,     4];
 
Problem.numconstraints = 2;
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = 1;
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = 1;

% 3. Call DIRECT
[fmin,x,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for GP test Function');
