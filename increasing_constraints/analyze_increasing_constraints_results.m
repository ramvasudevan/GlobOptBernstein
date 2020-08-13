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
% Updated: 7 Jan 2020
%
%% user parameters
% which problem to print
problem_index = 2 ; % pick an integer from 1 through 9

% colors
pcba_color = [0 0 1] ;
DIRECT_color = [140 140 40]./255 ;
bsos_color = [0 191 40]./255 ;
fmincon_color = [1 0 0] ;

% whether or not to plot everything
plot_figures_flag = true ;

% whether or not to save the output
save_pdfs_flag = false ;

%% data for reference (automated from here)
% names of problems
problem_name = {'ElAttar-Vidyasagar-Dutta','Powell','Wood',...
    'Dixon-Price 2-D','Dixon-Price 3-D','Dixon-Price 4-D',...
    'Beale','Bukin02','Deckkers-Aarts'} ;

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

% dimension and degree of problems
problem_dimension = [2 2 4 2 3 4 2 2 2] ;
problem_degree = [6 4 4 4 4 4 8 2 8] ;

% degree of constraints
constraint_degree = 2 ;

%% optimal value analysis (automated from here)
% total number of x ticks
total_steps = 20;

% number of constraints increase between x ticks
step = 10;

% load problem data
load(['more_',num2str(problem_index),'_info.mat']); % info about problem
load(['more_',num2str(problem_index),'.mat']); % problem cost, constraints, etc

% extract results and get error in optimal value
pcba_result = infos.bernstein_value_set - ground_truth(problem_index+1);
DIRECT_result = infos.DIRECT_value_set - ground_truth(problem_index+1);
bsos_result = infos.Lasserre_value_set - ground_truth(problem_index+1);
fmincon_result = infos.fmincon_value_set - ground_truth(problem_index+1);

% create labels
con_x = step:step:total_steps*step;

if plot_figures_flag
    % set up figure
    f1 = figure(1) ; clf ; hold on ;
    
    % plot fmincon boxes first
    boxplot_for_fmincon(fmincon_result,con_x,fmincon_color)
    
    % plot DIRECT result
    plot(1:total_steps,DIRECT_result,'o','Color',DIRECT_color,'MarkerSize',5,'LineWidth',3) ;
    
    % plot pcba result
    plot(1:total_steps,pcba_result,'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2) ;
    
    % plot bsos result
    plot(1:total_steps,bsos_result,'o','Color',bsos_color,'MarkerSize',10,'LineWidth',2) ;
    
    % labels
    xlabel('number of constraints \alpha','interpreter','tex');
    ylabel('(solver output) - (true optimum)','interpreter','tex');
    title(['''',problem_name{problem_index}, ''' solution error vs. number of constraints'],'interpreter','tex')
    xtickangle(45)
    set(gca,'FontSize',14)
    
    xtickangle(45)
    set(gca,'FontSize',14)
end

%% time analysis
% get data on time spent
pcba_time = infos.bernstein_time_set;
DIRECT_time = infos.DIRECT_time_set;
bsos_time = infos.Lasserre_time_set;
fmincon_time = infos.fmincon_time_set;

if plot_figures_flag
    % set up figure
    f2 = figure(2) ; clf ; hold on ;
    
    % plot fmincon box plots
    boxplot_for_fmincon(log10(fmincon_time),con_x,fmincon_color);
    
    % plot pcba result
    h_pcba = plot(1:total_steps,log10(pcba_time),'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2) ;
    
    % plot pcba result
    h_DIRECT = plot(1:total_steps,log10(DIRECT_time),'o','Color',DIRECT_color,'MarkerSize',5,'LineWidth',3) ;
    
    % plot bsos result
    h_bsos = plot(1:total_steps,log10(bsos_time),'o','Color',bsos_color,'MarkerSize',10,'LineWidth',2) ;
    
    % set axis bounds
    upp = max(log10(bsos_time));
    loww = log10(min(min(min(fmincon_time)),min(pcba_time)));
    axis([0 total_steps+1 loww upp+0.5]);
    
    % labels
    xlabel('number of constraints \alpha','interpreter','tex');
    ylabel('solve time [log_{10}(s)]','interpreter','tex');
    title(['''',problem_name{problem_index}, ''' solve time vs. number of constraints']) ;
    xtickangle(45)
    set(gca,'FontSize',14)
    
    % legend
    h_fmincon = findall(f2,'tag','Box') ;
    legend([h_pcba,h_DIRECT,h_bsos,h_fmincon(1)],{'PCBA','DIRECT','BSOS','fmincon'},'Location','northwest')
end

%% memory analysis
% get number of patches per number of constraints
number_of_items_per_number_of_constraints = infos.bernstein_apex_mem_set ;

% get memory per item in the list
cost = more_problem.cost' ;
memory_per_item = zeros(1,total_steps) ;
for idx = 1:total_steps
    % get constraints
    cons_cell = more_problem.constraints(1:10*idx) ;
    cons = cell2mat(cons_cell)' ;
    N_cons = 10*idx ;
    
    % get memory per item
    memory_per_item(idx) = get_memory_per_item(cost,cons,N_cons) ;
end

% peak memory usage
peak_memory_usage = number_of_items_per_number_of_constraints(:) .* memory_per_item(:) ./ 1024 ; % [kB]
mem_type = 'kB' ;
if max(peak_memory_usage) > 1024
    peak_memory_usage = peak_memory_usage ./ 1024 ;
    mem_type = 'MB' ;
end

if plot_figures_flag
    % set up figure
    f3 = figure(3) ; clf ;
    
    % plot pcba number of patches
    plot(con_x,number_of_items_per_number_of_constraints,'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2);
    ylabel('number of items');
    
    % labels
    xlabel('number of constraints \alpha','interpreter','tex');
    title(['''',problem_name{problem_index}, ''' number of items vs. number of constraints']) ;
    xticks(con_x)
    xtickangle(45)
    grid on ;
    set(gca,'FontSize',14)
    
    % set up figure for memory usage
    f4 = figure(4) ; clf ;
    
    % plot peak memory usage
    plot(con_x,peak_memory_usage,'x','Color',pcba_color,'MarkerSize',9,'LineWidth',2);
    
    % labels
    title(['''',problem_name{problem_index}, ''' peak memory usage vs. number of constraints']) ;
    ylabel(['peak memory usage [',mem_type,']'])
    xlabel('number of constraints \alpha','interpreter','tex');
    xticks(con_x)
    xtickangle(45)
    grid on ;
    set(gca,'FontSize',14)
end

%% display problem info
disp(' ')
disp(['--- ',upper(problem_name{problem_index}),' (INCREASING CONSTRAINTS) PCBA STATS ---'])
disp(['Decision variable dimension: ',num2str(size(cost,1)-1)]) ;
disp(['Max PCBA solve time spent: ',num2str(max(infos.bernstein_time_set))])
disp(['Max number of patches: ',num2str(max(infos.bernstein_apex_mem_set))])
disp(['Max GPU memory used: ',num2str(max(peak_memory_usage)),' ',mem_type])
disp(' ')

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
    
    % size the fourth figure correctly
    set(f4,'Units','Inches');
    pos = get(f4,'Position');
    set(f4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    % print
    print(f1,'increasing_constraints_error.pdf','-dpdf','-r0')
    print(f2,'increasing_constraints_time.pdf','-dpdf','-r0')
    print(f3,'increasing_constraints_pcba_number_of_items.pdf','-dpdf','-r0')
    print(f4,'increasing_constraints_pcba_memory_usage.pdf','-dpdf','-r0')
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
