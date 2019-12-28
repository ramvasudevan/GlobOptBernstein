clear; clc;
problem_index = 8;
%load(strcat('Benchmark_Evaluation/P',num2str(problem_index),'_infos.mat'),'infos');
eval(strcat('[raw_cost,raw_constraints] = setup_problem_matrix_P',num2str(problem_index),'();'));

%% Bernstein Algorithm
[bernstein_cost,bernstein_constraint,cons_length] = setup_problem_bernstein(raw_cost,raw_constraints);
[bernstein_opt,bernstein_memory,bernstein_accuracy] = PCBA(bernstein_cost,bernstein_constraint,cons_length,0,0);
% bernstein_memory = infos.bernstein_mem;
numIter = (min(find(bernstein_memory == 0,1,'first')) - 1) / 2;
plot(0.5:0.5:numIter,bernstein_memory(1:(numIter*2)),'.');
xlabel('iteration');
ylabel('number of patches');
title(['P',num2str(problem_index)]);

%% save output
save_filename = ['P',num2str(problem_index),'_time_and_memory_info.mat'] ;
save(save_filename,'numIter','bernstein_cost','bernstein_constraint','cons_length',...
    'bernstein_opt','bernstein_memory','bernstein_accuracy')

%% P7
% eval(strcat('[raw_cost,raw_constraints] = setup_problem_matrix_P7();'));
% P7_equalites = [4           0           0        0        -625;
%                 4           4           0        0      390625;
%                 0           4           1        0       -3125];
% numDimension = size(raw_cost,2) - 1;
% 
% %% Bernstein Algorithm
% [bernstein_cost,bernstein_constraint,cons_length] = setup_problem_bernstein(raw_cost,raw_constraints);
% [bernstein_opt,bernstein_memory,bernstein_accuracy] = PCBA(bernstein_cost,bernstein_constraint,cons_length,P7_equalites',[3]);
% numIter = (min(find(bernstein_memory == 0)) - 1) / 2;
% plot(bernstein_memory(1:(numIter*2)),'.');