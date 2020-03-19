% This is the script used to generate the randomized constraints in our
% paper V.C. All the plots included in this section of our paper 
% are based on the mat files stored in 'problem_matrix/', 
% named with 'more_X.mat'. If you would like to try a different problem,
% try to generate the problems in a different path.
%
% Author: Bohao Zhang
% Created: 6 Jan 2020

clear; clc;
total_steps = 20;
step = 10;

% change this flag to be true if you would like to try a different problem.
% Note that the original problem would be overwritten then 
% if you do not change the default file name .
save_problem_flag = false;
    
for problem_index = 0:9
    eval(['[cost,constraints,feasible_point] = setup_problem_matrix_more_',num2str(problem_index),'(',num2str(total_steps*step),');']);
    more_problem.cost = cost;
    more_problem.constraints = constraints;
    more_problem.feasible_point = feasible_point;
    
    if save_problem_flag == true
        save(['../problem_matrix/more_',num2str(problem_index),'.mat'],'more_problem');
    end
    
    clear more_problem;
end