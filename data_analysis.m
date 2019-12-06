%% optimal value analysis
% ba_filename = 'result/example_bernstein_result.mat';
% fmincon_filename = 'result/example_fmincon_result.mat';
% Lasserre_filename = 'result/example_Lasserre_result.mat';
problem_index = 7;
ba_filename = ['result/more_', num2str(problem_index),'_bernstein_result.mat'];
fmincon_filename = ['result/more_', num2str(problem_index),'_fmincon_result.mat'];
Lasserre_filename = ['result/more_', num2str(problem_index),'_Lasserre_result.mat'];
ba_result = load(ba_filename);
ba_result = ba_result.bernstein_value_set;
fmincon_result = load(fmincon_filename);
fmincon_result = fmincon_result.fmincon_value_set;
Lasserre_result = load(Lasserre_filename);
Lasserre_result = Lasserre_result.Lasserre_value_set;
total_steps = 20;
step = 10;
con_x = step:step:total_steps*step;

figure
f_hd = boxplot(fmincon_result,'Labels',con_x);
legend('fmincon')
hold on;
ba_hd = plot(1:total_steps,ba_result,'g*');
hold on;
bsos_hd = plot(1:total_steps,Lasserre_result,'mo');
xlabel('number of constraints');
ylabel('optimal value found');
legend([ba_hd, bsos_hd, f_hd(1,1)],'BA','BSOS','fmincon'); 
title('optimal value vs. number of constraints');

%% time analysis
% ba_filename = 'time/example_bernstein_time.mat';
% fmincon_filename = 'time/example_fmincon_time.mat';
% Lasserre_filename = 'time/example_Lasserre_time.mat';
ba_filename = ['time/more_', num2str(problem_index),'_bernstein_time.mat'];
fmincon_filename = ['time/more_', num2str(problem_index),'_fmincon_time.mat'];
Lasserre_filename = ['time/more_', num2str(problem_index),'_Lasserre_time.mat'];
ba_time = load(ba_filename);
ba_time = ba_time.bernstein_time_set;
fmincon_time = load(fmincon_filename);
fmincon_time = fmincon_time.fmincon_time_set;
Lasserre_time = load(Lasserre_filename);
Lasserre_time = Lasserre_time.Lasserre_time_set;

figure
f_hd = boxplot(log(fmincon_time),'Labels',con_x);
hold on;
ba_hd = plot(1:total_steps,log(ba_time),'g*');
hold on;
bsos_hd = plot(1:total_steps,log(Lasserre_time),'mo');
upp = max(log(Lasserre_time));
loww = log(min(min(min(fmincon_time)),min(ba_time)));
axis([0 total_steps+1 loww upp+0.5]);
xlabel('number of constraints');
ylabel('solve time log(s)');
legend([ba_hd, bsos_hd, f_hd(1,1)],'BA','BSOS','fmincon'); 
title('time vs. number of constraints');