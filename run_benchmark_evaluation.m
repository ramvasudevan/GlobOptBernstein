%% Introduction
% This is the script that runs the results for benchmark problems P1-P8 (except P7) in
% our paper V.B. All the results all stored in 'Benchmark_Evaluation/'
% named with 'PX_infos.mat' and they are consistent with what is shown 
% in this section of our paper. If you want a different comparison, 
% try to save the results in a different path.

%% setup the problem
clear all; clc;
problem_index = 8;
% the true optimal value (or the best we found) for each problem
ground_truth = [-5.5080132636,...
                -6961.8138816446,...
                3.0000011115,...
                -4,...
                0,...
                6395.5078,...
                1.0898639714,...
                42.4440570797];
eval(strcat('[raw_cost,raw_constraints] = setup_problem_matrix_P',num2str(problem_index),'();'));
numDimension = size(raw_cost,2) - 1;

%% Bernstein Algorithm
[bernstein_cost,bernstein_constraint,cons_length] = setup_problem_bernstein(raw_cost,raw_constraints);

% run for memory analysis
verboseMode = 0;
memoryRecordMode = 1;
pcba_options = memoryRecordMode * 2 + verboseMode;
[bernstein_opt,bernstein_accuracy,bernstein_memory] = PCBA(bernstein_cost,bernstein_constraint,cons_length,0,0,pcba_options);

% run another 50 times to get an average time
pcba_num = 50;
verboseMode = 0;
memoryRecordMode = 0;
pcba_options = memoryRecordMode * 2 + verboseMode;
bernstein_time_set = nan(pcba_num,1);
for num = 1:pcba_num
    bernstein_start_t = tic;
    [bernstein_opt,bernstein_accuracy] = PCBA(bernstein_cost,bernstein_constraint,cons_length,0,0,pcba_options);
    bernstein_time_set(num) = toc(bernstein_start_t);
end

if bernstein_opt == -12345
    bernstein_exitflag = -1;
else
    bernstein_opt = bernstein_opt(:,1);
    bernstein_exitflag = 1;
    [bernstein_value,bernstein_feasibility,bernstein_violate_terms,bernstein_difference] = evaluate_opt_result(raw_cost,raw_constraints,bernstein_opt);
end

%% fmincon
clc;
fmincon_num = 50;
fmincon_time_set = nan(fmincon_num,1);
fmincon_value_set = nan(fmincon_num,1);
fmincon_exit_flag = nan(fmincon_num, 1); % 1 means successful, 0 means not
fmincon_cost = @(k) setup_cost_fmincon(raw_cost,k);
fmincon_nonlcon = @(k) setup_constraints_fmincon(raw_constraints,{},k);

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

%% Lasserre
Lasserre_number = 50;
BSOSsolver = 'sqlp';
SBSOSsolver = 'sqlp';
Lasserre_d_choice = [3,2,2,4,2,3,0,2];
Lasserre_k_choice = [3,2,2,4,2,3,0,2];
Lasserre_d = Lasserre_d_choice(problem_index);
Lasserre_k = Lasserre_k_choice(problem_index);
Lasserre_time_set = nan(Lasserre_number,1);
Lss_constraints = scale_for_lss(raw_constraints);
setup_Lasserre = @()setup_problem_Lasserre(raw_cost,Lss_constraints,Lasserre_d,Lasserre_k);
tag = 'setup_Lasserre';
for i = 1:Lasserre_number
    clc;
    Lasserre_start_t = tic;
    eval(['[pop.F,pop.G,pop.I,pop.J,pop.d,pop.k] = ',tag,'();']);
    k = pop.k;
    pop.n = size(pop.F,2)-1;
    psol_temp = NaN;
    for d = 1:pop.d
        %% BSOS
        algo = 'BSOS';
        solver = eval([algo,'solver',';']);
        pop.d = d; pop.k = k;
        psol_temp = lss(pop,tag,algo,solver);
        clear sdp sol psol;

        %% SBSOS
    %         algo = 'SBSOS';
    %         solver = eval([algo,'solver',';']);
    %         pop.d = d; pop.k = k;
    %         psol_temp = lss(pop,tag,algo,solver);
    %         clear sdp sol psol;
    end
    Lasserre_time_set(i) = toc(Lasserre_start_t);
    Lasserre_value = psol_temp.obj;
end

%% DIRECT
DIRECT_number = 50;
DIRECT_time_set = nan(DIRECT_number,1);
options.maxevals  = 5000;
options.maxits    = 100;
options.testflag  = 0;
options.showits   = 1;
options.tol       = bernstein_accuracy;

for i = 1:DIRECT_number
    clc;
    DIRECT_start_t = tic;
    eval(strcat('[DIRECT_result,DIRECT_opt] = DIRECT_Benchmark_P',num2str(problem_index),'(options);'));
    DIRECT_time_set(i) = toc(DIRECT_start_t);
    [DIRECT_value,DIRECT_feasibility,DIRECT_violate_terms,DIRECT_difference] = evaluate_opt_result(raw_cost,raw_constraints,DIRECT_opt);
end

%% save the data
infos.fmincon_time_set = fmincon_time_set;
infos.fmincon_value_set = fmincon_value_set;
infos.fmincon_exit_flag = fmincon_exit_flag;
infos.bernstein_time_set = bernstein_time_set;
infos.bernstein_value = bernstein_value;
infos.bernstein_mem = bernstein_memory;
infos.bernstein_accuracy = bernstein_accuracy;
infos.Lasserre_time_set = Lasserre_time_set;
infos.Lasserre_value = Lasserre_value;
infos.Lasserre_d = Lasserre_d;
infos.Lasserre_k = Lasserre_k;
infos.DIRECT_time_set = DIRECT_time_set;
infos.DIRECT_value = DIRECT_value;
save(strcat('Benchmark_Evaluation/P',num2str(problem_index),'_infos.mat'),'infos');

%% data analysis
disp('PCBA time median:')
disp(median(infos.bernstein_time_set));
disp('PCBA error:')
disp(infos.bernstein_value - ground_truth(problem_index));
disp('PCBA stopping crtieria:')
disp(infos.bernstein_accuracy);
disp(' ')
disp('BSOS time median:')
disp(median(infos.Lasserre_time_set));
disp('BSOS error:')
disp(infos.Lasserre_value - ground_truth(problem_index));
disp('BSOS d & k choice:')
disp([infos.Lasserre_d,infos.Lasserre_k]);
disp(' ')
disp('fmincon time median:')
disp(median(infos.fmincon_time_set));
disp('fmincon error median:')
disp(median(infos.fmincon_value_set) - ground_truth(problem_index));
disp('DIRECT time median:')
disp(median(infos.DIRECT_time_set));
disp('DIRECT error median:')
disp(median(infos.DIRECT_value) - ground_truth(problem_index));
