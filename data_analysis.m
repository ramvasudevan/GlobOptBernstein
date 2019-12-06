%% optimal value analysis
problem_index = 0;
ground_truth = [0;0;-124.75;-24777;0;0;0;1.712780354;0;0];
total_steps = 20;
step = 10;
load(['Increasing_Number_of_Constraints/more_',num2str(problem_index),'_info.mat']);
ba_result       = infos.bernstein_value_set - ground_truth(problem_index+1);
fmincon_result  = infos.fmincon_value_set   - ground_truth(problem_index+1);
% Lasserre_result = infos.Lasserre_value_set  - ground_truth(problem_index+1);

con_x = step:step:total_steps*step;

figure
f_hd = boxplot(fmincon_result,'Labels',con_x);
legend('fmincon')
hold on;
ba_hd = plot(1:total_steps,ba_result,'g*');
hold on;
% bsos_hd = plot(1:total_steps,Lasserre_result,'mo');
xlabel('number of constraints');
ylabel('optimal value found');
title('optimal value vs. number of constraints');

%% time analysis
ba_time       = infos.bernstein_time_set;
fmincon_time  = infos.fmincon_time_set;
% Lasserre_time = infos.Lasserre_time_set;

figure
f_hd = boxplot(log(fmincon_time),'Labels',con_x);
hold on;
ba_hd = plot(1:total_steps,log(ba_time),'g*');
hold on;
% bsos_hd = plot(1:total_steps,log(Lasserre_time),'mo');
upp = max(log(Lasserre_time));
loww = log(min(min(min(fmincon_time)),min(ba_time)));
axis([0 total_steps+1 loww upp+0.5]);
xlabel('number of constraints');
ylabel('solve time log(s)');
title('time vs. number of constraints');