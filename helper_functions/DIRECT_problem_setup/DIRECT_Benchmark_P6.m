function [fmin,x] = DIRECT_Benchmark_P6(options)
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
    options.maxits    = 30;
    options.testflag  = 0;
    options.showits   = 1;
    options.tol = 4.2677e-04 ;
end

% 1. Establish bounds for variables
bounds = [0 1;0 1;0 1;0 1];
constraint_penalty = 1 ;

% 2. Send options to Direct

% 2a. NEW!
% Pass Function as part of a Matlab Structure
cost = [[   0,   0,   0,   0, 6395.507828125]
        [ 1.0,   0,   0,   0, 567.11175]
        [ 2.0,   0,   0,   0, 40.070953125]
        [   0, 1.0,   0,   0, 1504.439296875]
        [   0,   0, 1.0,   0, 907.1534375]
        [   0,   0, 2.0,   0, 27.7828125]
        [ 1.0,   0, 1.0,   0, 37.2]
        [   0, 1.0, 1.0,   0, 316.7240625]
        [   0, 1.0, 2.0,   0, 16.6696875]
        [   0,   0,   0, 1.0, 720.0622]
        [ 1.0,   0,   0, 1.0, 52.24065]
        [ 2.0,   0,   0, 1.0, 9.795121875]
        [   0,   0, 1.0, 1.0, 68.464]];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 4;

g1 = [[     0,   0,   0,   0, -0.08325]
        [ 1.0,   0,   0,   0, -0.375]
        [   0,   0,   1.0, 0,  0.0965]];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = constraint_penalty;
  
g2 = [[   0,   0,   0,   0, -0.17185]
    [     0, 1.0,   0,   0,  -0.375]
    [     0,   0, 1.0,   0,  0.0477]];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = constraint_penalty;

g3 = [[   0,   0,   0,     0, -1086109.9856481687165796756744385]
        [ 0,   0,   1.0,   0, -276067.45443420308082990478730569]
        [ 0,   0,   2.0,   0, -21991.148575128552669238503682957]
        [ 0,   0,   3.0,   0, -523.59877559829887307710723054658]
        [ 0,   0,   0,   1.0, -155940.80534256336187418946093754]
        [ 0,   0,   1.0, 1.0, -32829.643230013339341934623355271]
        [ 0,   0,   2.0, 1.0, -1727.8759594743862811544538608037]];
Problem.constraint(3).func = @(x)evaluate_function(g3,x);
Problem.constraint(3).penalty = constraint_penalty;

g4 = [0  0  0  0  -150;
      0  0  0  1    22];
Problem.constraint(4).func = @(x)evaluate_function(g4,x);
Problem.constraint(4).penalty = constraint_penalty;

% 3. Call DIRECT
[fmin,x,~] = Direct(Problem,bounds,options);

end

