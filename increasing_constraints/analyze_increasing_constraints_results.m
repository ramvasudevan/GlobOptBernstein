%% description
% This script plots the data for the "increasing constraints" problems that
% compare our PCBA algorithm, Lasserre's BSOS, and MATLAB's fmincon.
%
% The problem_index variable determines which problem to plot:
%   1 - El-Attar-Vidyasagar-Dutta
%   2 - Powell
%   3 - Wood
%   4 - Dixon-Price 2-D
%   5 - Dixon-Price 3-D
%   6 - Dixon-Price 4-D
%   7 - Beale
%   8 - Bukin02
%   9 - Deckkers-Aarts
%
% Authors: Bohao Zhang and Shreyas Kousik
% Created: Summer 2019
% Updated: 4 Jan 2020
%
%% user parameters
% which problem to print
problem_index = 8 ; % pick an integer from 1 through 9

% colors
pcba_color = [0 0 1] ;
bsos_color = [0 191 40]./255 ;
fmincon_color = [1 0 0] ;

% whether or not to save the output
save_pdfs_flag = false ;

%% optimal value analysis (automated from here)
% true optimal value for each cost function
ground_truth = [0;
                0;
                -124.75;
                -24.77109375;
                0;
                0;
                0;
                0.000001712780354;
                0;
                0];

% total number of x ticks
total_steps = 20;

% number of constraints increase between x ticks
step = 10;

% load data
load(['more_',num2str(problem_index),'_info.mat']);

% extract results and get error in optimal value
pcba_result = infos.bernstein_value_set - ground_truth(problem_index+1);
bsos_result = infos.Lasserre_value_set - ground_truth(problem_index+1);
fmincon_result = infos.fmincon_value_set - ground_truth(problem_index+1);

% create labels
con_x = step:step:total_steps*step;

% set up figure
f1 = figure(1) ; clf ; hold on ;

% plot fmincon boxes first
boxplot_for_fmincon(fmincon_result,con_x,fmincon_color)

% plot pcba result
plot(1:total_steps,pcba_result,'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2) ;

% plot bsos result
plot(1:total_steps,bsos_result,'o','Color',bsos_color,'MarkerSize',10,'LineWidth',2) ;

% labels
xlabel('number of constraints');
ylabel('(solver output) - (true optimum)');
title('solution error vs. number of constraints');
xtickangle(45)
set(gca,'FontSize',14)

xtickangle(45)
set(gca,'FontSize',14)

%% time analysis
% get data on time spent
pcba_time = infos.bernstein_time_set;
bsos_time = infos.Lasserre_time_set;
fmincon_time = infos.fmincon_time_set;

% set up figure
f2 = figure(2) ; clf ; hold on ;

% plot fmincon box plots
boxplot_for_fmincon(log10(fmincon_time),con_x,fmincon_color);

% plot pcba result
h_pcba = plot(1:total_steps,log10(pcba_time),'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2) ;

% plot bsos result
h_bsos = plot(1:total_steps,log10(bsos_time),'o','Color',bsos_color,'MarkerSize',10,'LineWidth',2) ;

% set axis bounds
upp = max(log10(bsos_time));
loww = log10(min(min(min(fmincon_time)),min(pcba_time)));
axis([0 total_steps+1 loww upp+0.5]);

% labels
xlabel('number of constraints');
ylabel('solve time [log_{10}(s)]');
title('solve time vs. number of constraints') ;
xtickangle(45)
set(gca,'FontSize',14)

% legend
h_fmincon = findall(f2,'tag','Box') ;
legend([h_pcba,h_bsos,h_fmincon(1)],{'PCBA','BSOS','fmincon'},'Location','northwest')

%% memory analysis
bernstein_apex_mem = infos.bernstein_apex_mem_set;

% set up figure
f3 = figure(3) ; clf ; hold on ;

% plot pcba memory usage
plot(con_x,bernstein_apex_mem,'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2);

% labels
xlabel('number of constraints');
ylabel('number of patches]');
title('number of patches vs. number of constraints') ;
xtickangle(45)
set(gca,'FontSize',14)

%% save output
if save_pdfs_flag
    % size the first figure correctly
    set(f1,'Units','Inches');
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    % size the second figure correctly
    set(f2,'Units','Inches');
    pos = get(f2,'Position');
    set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    % size the third figure correctly
    set(f3,'Units','Inches');
    pos = get(f3,'Position');
    set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print
    print(f1,'increasing_constraints_error.pdf','-dpdf','-r0')
    print(f2,'increasing_constraints_time.pdf','-dpdf','-r0')
    print(f3,'increasing_constraints_pcba_memory.pdf','-dpdf','-r0')
end

%% helper function
function boxplot_for_fmincon(data,labels,color)
    % make the boxplot
    boxplot(data,'Labels',labels) ;

    % set the median and box lines to red
    m = findobj(gcf,'type','Line','Tag','Median') ;
    b = findobj(gcf,'type','Line','Tag','Box') ;
    set(m,'Color',color)
    set(b,'Color',color) ;

    % set the linewidths to 1
    l = findobj(gcf,'type','Line') ;
    set(l,'LineWidth',1)

    % set the markers to have thicker lines
    k = findobj(gcf,'tag','Outliers') ;
    set(k,'LineWidth',1.5)
end
