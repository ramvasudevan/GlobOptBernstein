function [fmin,x] = DIRECT_Benchmark_P8(options)
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
cost = [[   0      0      0      0   35.784]
        [ 1.0      0      0      0   57.936]
        [   0    1.0      0      0  177.216]
        [   0      0    1.0      0   42.600]
        [ 1.0      0    1.0      0  289.680]
        [   0      0      0    1.0  102.240]
        [   0    1.0      0    1.0  708.864]
        [   0      0    1.0    1.0  -34.080]];
Problem.f = @(x)evaluate_function(cost,x);

Problem.numconstraints = 7;

g1 = [[   0,   0,   0,   0,     6.39453125]
        [ 1.0,   0,   0,   0,  -107.046875]
        [ 2.0,   0,   0,   0,   -1029.5625]
        [ 3.0,   0,   0,   0,     -1228.25]
        [   0, 1.0,   0,   0,   -80.640625]
        [ 1.0, 1.0,   0,   0,    -953.0625]
        [ 2.0, 1.0,   0,   0,     -2817.75]
        [   0,   0, 1.0,   0, -49.62890625]
        [   0,   0, 2.0,   0,  22.55859375]
        [   0,   0, 3.0,   0,  -3.41796875]
        [ 1.0,   0, 1.0,   0,   -613.59375]
        [ 2.0,   0, 1.0,   0,   -1896.5625]
        [ 1.0,   0, 2.0,   0,   139.453125]
        [   0, 1.0, 1.0,   0,  -368.671875]
        [   0, 1.0, 2.0,   0,   167.578125]
        [   0, 1.0, 3.0,   0,   -25.390625]
        [ 1.0, 1.0, 1.0,   0,    -4558.125]
        [ 2.0, 1.0, 1.0,   0,    -14088.75]
        [ 1.0, 1.0, 2.0,   0,    1035.9375]
        [   0,   0,   0, 1.0,   -20.796875]
        [ 1.0,   0,   0, 1.0,    -385.6875]
        [ 2.0,   0,   0, 1.0,     -2384.25]
        [ 3.0,   0,   0, 1.0,      -4913.0]
        [   0,   0, 1.0, 1.0,    28.359375]
        [   0,   0, 2.0, 1.0,   -12.890625]
        [   0,   0, 3.0, 1.0,     1.953125]
        [ 1.0,   0, 1.0, 1.0,      350.625]
        [ 2.0,   0, 1.0, 1.0,      1083.75]
        [ 1.0,   0, 2.0, 1.0,     -79.6875]];
Problem.constraint(1).func = @(x)evaluate_function(g1,x);
Problem.constraint(1).penalty = 1;
  
g2 = [[   0,   0,   0,   0,                        43.41080711]
        [ 1.0,   0,   0,   0,                      -243.046875]
        [ 2.0,   0,   0,   0,                       -1029.5625]
        [ 3.0,   0,   0,   0,                         -1228.25]
        [   0, 1.0,   0,   0,                       -80.640625]
        [ 1.0, 1.0,   0,   0,                        -953.0625]
        [ 2.0, 1.0,   0,   0,                         -2817.75]
        [   0,   0, 1.0,   0,                     -49.62890625]
        [   0,   0, 2.0,   0,                      22.55859375]
        [   0,   0, 3.0,   0,                      -3.41796875]
        [ 1.0,   0, 1.0,   0,                       -613.59375]
        [ 2.0,   0, 1.0,   0,                       -1896.5625]
        [ 1.0,   0, 2.0,   0,                       139.453125]
        [   0, 1.0, 1.0,   0,                      -368.671875]
        [   0, 1.0, 2.0,   0,                       167.578125]
        [   0, 1.0, 3.0,   0,                       -25.390625]
        [ 1.0, 1.0, 1.0,   0,                        -4558.125]
        [ 2.0, 1.0, 1.0,   0,                        -14088.75]
        [ 1.0, 1.0, 2.0,   0,                        1035.9375]
        [   0,   0,   0, 1.0,                       -20.796875]
        [ 1.0,   0,   0, 1.0,                        -385.6875]
        [ 2.0,   0,   0, 1.0,                         -2384.25]
        [ 3.0,   0,   0, 1.0,                          -4913.0]
        [   0,   0, 1.0, 1.0,                        28.359375]
        [   0,   0, 2.0, 1.0,                       -12.890625]
        [   0,   0, 3.0, 1.0,                         1.953125]
        [ 1.0,   0, 1.0, 1.0,                          350.625]
        [ 2.0,   0, 1.0, 1.0,                          1083.75]
        [ 1.0,   0, 2.0, 1.0,                         -79.6875]];
Problem.constraint(2).func = @(x)evaluate_function(g2,x);
Problem.constraint(2).penalty = 1;

g3 = [[   0,   0,   0,   0,    -5.382080078125]
        [ 1.0,   0,   0,   0,   -86.1455078125]
        [ 2.0,   0,   0,   0,    -414.30859375]
        [ 3.0,   0,   0,   0,      -537.359375]
        [   0, 1.0,   0,   0,   -26.3427734375]
        [ 1.0, 1.0,   0,   0,    -361.71484375]
        [ 2.0, 1.0,   0,   0,     -1232.765625]
        [   0,   0, 1.0,   0,    -55.498046875]
        [   0,   0, 2.0,   0,  -98.69384765625]
        [   0,   0, 3.0,   0,       47.8515625]
        [   0,   0, 4.0,   0,  -7.476806640625]
        [ 1.0,   0, 1.0,   0,  -741.6748046875]
        [ 2.0,   0, 1.0,   0,    -2901.2890625]
        [ 3.0,   0, 1.0,   0,     -2686.796875]
        [ 1.0,   0, 2.0,   0, -1281.2255859375]
        [ 2.0,   0, 2.0,   0,   -4148.73046875]
        [ 1.0,   0, 3.0,   0,   305.0537109375]
        [   0, 1.0, 1.0,   0,     -337.6953125]
        [   0, 1.0, 2.0,   0,   -733.154296875]
        [   0, 1.0, 3.0,   0,        355.46875]
        [   0, 1.0, 4.0,   0,   -55.5419921875]
        [ 1.0, 1.0, 1.0,   0,   -4079.00390625]
        [ 2.0, 1.0, 1.0,   0,     -12327.65625]
        [ 1.0, 1.0, 2.0,   0,   -9517.67578125]
        [ 2.0, 1.0, 2.0,   0,    -30819.140625]
        [ 1.0, 1.0, 3.0,   0,    2266.11328125]
        [   0,   0,   0, 1.0,    -4.7861328125]
        [   0,   0,   0, 2.0,           -1.875]
        [ 1.0,   0,   0, 1.0,    -138.98828125]
        [ 2.0,   0,   0, 1.0,     -1043.109375]
        [ 3.0,   0,   0, 1.0,       -2149.4375]
        [   0, 1.0,   0, 1.0,             32.5]
        [   0, 1.0,   0, 2.0,            -13.0]
        [ 1.0, 1.0,   0, 1.0,            221.0]
        [   0,   0, 1.0, 1.0,      -36.5234375]
        [   0,   0, 2.0, 1.0,     56.396484375]
        [   0,   0, 3.0, 1.0,        -27.34375]
        [   0,   0, 4.0, 1.0,     4.2724609375]
        [   0,   0, 1.0, 2.0,            0.625]
        [ 1.0,   0, 1.0, 1.0,    -711.54296875]
        [ 2.0,   0, 1.0, 1.0,      -4741.40625]
        [ 3.0,   0, 1.0, 1.0,      -10747.1875]
        [ 1.0,   0, 2.0, 1.0,     732.12890625]
        [ 2.0,   0, 2.0, 1.0,      2370.703125]
        [ 1.0,   0, 3.0, 1.0,    -174.31640625]];
Problem.constraint(3).func = @(x)evaluate_function(g3,x);
Problem.constraint(3).penalty = 1;

g4 = [   0     0    0    0   -3;
         1     0    0    0    17;
         0     1    0    0   -39];
Problem.constraint(4).func = @(x)evaluate_function(g4,x);
Problem.constraint(4).penalty = 1;

g5 = [0     0     0    0    1;
      1     0     0    0   -17;
      0     1     0    0    26];
Problem.constraint(5).func = @(x)evaluate_function(g5,x);
Problem.constraint(5).penalty = 1;

g6 = [0     0     0    0   -0.2500;
      0     0     1    0    0.6250;
      0     0     0    1   -1.5000];
Problem.constraint(6).func = @(x)evaluate_function(g6,x);
Problem.constraint(6).penalty = 1;

g7 = [0     0     1    0   -0.6250
      0     0     0    1    0.5000];
Problem.constraint(7).func = @(x)evaluate_function(g7,x);
Problem.constraint(7).penalty = 1;

% 3. Call DIRECT
[fmin,x,~] = Direct(Problem,bounds,options);

end

