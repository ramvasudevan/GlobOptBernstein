%% Introduction
% This is the script that runs the results for benchmark problems P7 in
% our paper V.B. All the results all stored in 'Benchmark_Evaluation/'
% named with 'P7_infos.mat' and they are consistent with what is shown 
% in this section of our paper. If you want a different comparison, 
% try to save the results in a different path.

%% setup the problem
clc;
ground_truth = 1.089;
fmincon_num = 50;
fmincon_time_set = nan(fmincon_num, 1);
fmincon_value_set = nan(fmincon_num, 1);
fmincon_exit_flag = nan(fmincon_num, 1); % 1 means successful, 0 means not

eval(strcat('[raw_cost,raw_constraints] = setup_problem_matrix_P7();'));
P7_equalites = [4           0           0        0        -625;
                4           4           0        0      390625;
                0           4           1        0       -3125];
numDimension = size(raw_cost,2) - 1;

%% Bernstein Algorithm
[bernstein_cost,bernstein_constraint,cons_length] = setup_problem_bernstein(raw_cost,raw_constraints);
bernstein_start_t = tic;
[bernstein_opt,bernstein_apex_memory,bernstein_accuracy] = bernstein(bernstein_cost,bernstein_constraint,cons_length,P7_equalites',[3]);
bernstein_time = toc(bernstein_start_t);
if bernstein_opt == -12345
    bernstein_exitflag = -1;
else
    bernstein_opt = bernstein_opt(:,1);
    bernstein_exitflag = 1;
    [bernstein_value,bernstein_feasibility,bernstein_violate_terms,bernstein_difference] = evaluate_opt_result(raw_cost,raw_constraints,bernstein_opt);
end

%% fmincon
fmincon_cost = @(k) setup_cost_fmincon(raw_cost,k);
fmincon_nonlcon = @(k) setup_constraints_fmincon(raw_constraints,{P7_equalites},k);

% create optimization options
options =  optimoptions('fmincon',...
                'MaxFunctionEvaluations',1e5,...
                'MaxIterations',1e4,...
                'OptimalityTolerance',bernstein_accuracy,...
                'CheckGradients',false,...
                'FiniteDifferenceType','central',...
                'Diagnostics','off',...
                'SpecifyConstraintGradient',true,...
                'SpecifyObjectiveGradient',true);

% call fmincon
for j = 1:fmincon_num
    fmincon_start_t = tic;
    initial_guess = rand(size(raw_cost,2)-1,1);
    try
        [fmincon_opt,~,fmincon_exitflag] = fmincon(fmincon_cost,...
                                    initial_guess,...
                                    [],[],... % linear inequality constraints
                                    [],[],... % linear equality constraints
                                    zeros(1,numDimension),... % lower bounds
                                    ones(1,numDimension),... % upper bounds
                                    fmincon_nonlcon,...
                                    options) ;
    catch
        fmincon_exitflag = -1 ;
    end
    fmincon_time = toc(fmincon_start_t);
    fmincon_time_set(j) = fmincon_time;
    if fmincon_exitflag ~= -1
        [fmincon_value,fmincon_feasibility,fmincon_violate_terms,fmincon_difference] = evaluate_opt_result(raw_cost,raw_constraints,fmincon_opt);
        fmincon_value_set(j) = fmincon_value;
        fmincon_exit_flag(j) = 1;
    else
        fmincon_exit_flag(j) = 0;
    end
end

%% save the data
infos.fmincon_time_set = fmincon_time_set;
infos.fmincon_value_set = fmincon_value_set;
infos.fmincon_exit_flag = fmincon_exit_flag;
infos.bernstein_time = bernstein_time;
infos.bernstein_value = bernstein_value;
infos.bernstein_apex_mem = bernstein_apex_memory;
% save(strcat('Benchmark_Evaluation/P7_infos.mat'),'infos');

%% data analysis
disp(median(infos.fmincon_time_set));
disp(median(infos.fmincon_value_set) - ground_truth);
disp(infos.bernstein_time);
disp(infos.bernstein_value - ground_truth);