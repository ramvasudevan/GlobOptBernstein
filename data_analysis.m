%% optimal value analysis
problem_index = 2 ; % 2 is good
ground_truth = [0;0;-124.75;-24.77109375;0;0;0;0.000001712780354;0;0];
total_steps = 20;
step = 10;
load(['Increasing_Number_of_Constraints/more_',num2str(problem_index),'_info.mat']);
ba_result       = infos.bernstein_value_set - ground_truth(problem_index+1);
fmincon_result  = infos.fmincon_value_set   - ground_truth(problem_index+1);
Lasserre_result = infos.Lasserre_value_set  - ground_truth(problem_index+1);

con_x = step:step:total_steps*step;

figure(1) ; clf ;
f_hd = boxplot(fmincon_result,'Labels',con_x);
hold on;
ba_hd = plot(1:total_steps,ba_result,'bx','MarkerSize',9,'LineWidth',2);
hold on;
bsos_hd = plot(1:total_steps,Lasserre_result,'ro','MarkerSize',10,'LineWidth',1.5);
xlabel('number of constraints');
ylabel('error');
title('error vs. number of constraints');

xtickangle(45)
set(gca,'FontSize',14)

%% time analysis
ba_time       = infos.bernstein_time_set;
fmincon_time  = infos.fmincon_time_set;
Lasserre_time = infos.Lasserre_time_set;

figure(2) ; clf ;
f_hd = boxplot(log(fmincon_time),'Labels',con_x);
hold on;
ba_hd = plot(1:total_steps,log(ba_time),'bx','MarkerSize',9,'LineWidth',2);
hold on;
bsos_hd = plot(1:total_steps,log(Lasserre_time),'ro','MarkerSize',10,'LineWidth',1.5);
upp = max(log(Lasserre_time));
loww = log(min(min(min(fmincon_time)),min(ba_time)));
axis([0 total_steps+1 loww upp+0.5]);
xlabel('number of constraints');
ylabel('solve time log(s)');
title('time vs. number of constraints') ;

xtickangle(45)
set(gca,'FontSize',14)