clear; clc;

%% 1. Go to src/ and compile the mexcuda program

%% 2. Set up the problem
% Problem: min cost, s.t.
%              g_i >= 0
%
% Express the polynomial using the following matrix form
% [[degree of x1, degree of x2, ..., degree of xn, coefficient], (the first monomial in the polynomial)
%  [degree of x1, degree of x2, ..., degree of xn, coefficient], (the second monomial in the polynomial)
%  ...
%  [degree of x1, degree of x2, ..., degree of xn, coefficient]] (the last monomial in the polynomial)

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

g1 = [[     0,   0,   0,   0, -0.08325]
        [ 1.0,   0,   0,   0, -0.375]
        [   0,   0,   1.0, 0,  0.0965]];
  
g2 = [[   0,   0,   0,   0, -0.17185]
    [     0, 1.0,   0,   0,  -0.375]
    [     0,   0, 1.0,   0,  0.0477]];
       
g3 = [[   0,   0,   0,     0, -1086109.9856481687165796756744385]
        [ 0,   0,   1.0,   0, -276067.45443420308082990478730569]
        [ 0,   0,   2.0,   0, -21991.148575128552669238503682957]
        [ 0,   0,   3.0,   0, -523.59877559829887307710723054658]
        [ 0,   0,   0,   1.0, -155940.80534256336187418946093754]
        [ 0,   0,   1.0, 1.0, -32829.643230013339341934623355271]
        [ 0,   0,   2.0, 1.0, -1727.8759594743862811544538608037]];

g4 = [0  0  0  0  -150;
      0  0  0  1    22];
  
constraints = {g1;g2;g3;g4};

% convert constraints cell to a larger matrix
[bernstein_cost,bernstein_constraint,cons_length] = setup_problem_bernstein(cost,constraints);

%% 3. Run the optimization program
% Program options
verboseMode = 0; % turn on/off the iteration output
memoryRecordMode = 1; % record the memory usage (patch number) or not
pcba_options = memoryRecordMode * 2 + verboseMode;

% go to src/BCDATA.h for further internal parameters infomation. You will
% have to recompile the program if you want to change them

% PCBA
% Input: 
% -- cost matrix
% -- total constraint matrix
% -- number of rows of each constraint matrix in the total constraint matrix
% -- total equality constraint matrix
% -- number of rows of each equality constraint matrix in the total total equality constraint matrix
% -- output options
% Output:
% -- optimization result
% -- accuracy achieved (approximated difference to actual optimum)
% -- memory usage at each subdivision iteration
[bernstein_opt,bernstein_accuracy,bernstein_memory] = PCBA(bernstein_cost,bernstein_constraint,cons_length,0,0,pcba_options);

%% 4. Evaluate the result
% evaluate_opt_result
% Input: 
% -- cost matrix
% -- constraint matrix cell
% -- optimization result
% Output:
% -- result optimum value
% -- result optimum feasibility
% -- index of constraints that the result optimum violates
% -- difference in constraints that the result optimum violates
[bernstein_value,bernstein_feasibility,bernstein_violate_terms,bernstein_difference] = evaluate_opt_result(cost,constraints,bernstein_opt);

