%% description
% This script plots the number of Bernstein patches (i.e., the length of
% the list) versus the number of iterations of PCBA for each of the 8
% benchmark problems.
%
% Author: Shreyas Kousik and Bohao Zhang
% Created: 28 Dec 2019
% Updated: 7 Jan 2020
%
%% user parameters
% which problem to plot
problem_index = 4 ;

% whether or not to save the output
save_pdf_flag = false ;

% the maximum number of iteration used in the program
default_max_iteration = 28;

%% automated from here
% load data
load(['P',num2str(problem_index),'_infos.mat'])

%% process data
% get the problem dimension and degree in each dimension ;
eval(strcat('[raw_cost,raw_constraints] = setup_problem_matrix_P',num2str(problem_index),'();'));
[bernstein_cost,bernstein_constraint,cons_length] = setup_problem_bernstein(raw_cost,raw_constraints);
dimension = size(bernstein_cost,1) - 1 ;
degrees = max([bernstein_cost(1:dimension,:),bernstein_constraint(1:dimension,:)],[],2) ;

% relabel the "bernstein memory" variable
bernstein_N_patches = infos.bernstein_mem ;

% (over)approximate the memory used
N_cons = length(cons_length) ; % N_cons + 1 is # of polynomials represented by items in list
memory_per_item = get_memory_per_item(bernstein_cost,bernstein_constraint,N_cons) ;
bernstein_memory = memory_per_item.*bernstein_N_patches ; % in bytes

% get mem usage in kB or mB
mem_type = 'kB' ;
memory_per_item = memory_per_item ./ 1024 ;
bernstein_memory = bernstein_memory ./ 1024 ;

if max(bernstein_memory) > 1000
    mem_type = 'MB' ;
    memory_per_item = memory_per_item ./ 1024 ;
    bernstein_memory = bernstein_memory ./ 1024 ;
end

% get the number of iterations (each entry in bernstein_memory is the
% number of patches either from subdivision or elimination, and subdivision
% happens once along each dimension in each iteration, so the denominator
% is 2*d; the first "0" entry in bernstein_memory marks the iteration
% where the problem was solved to the specified tolerances)
problem_solved_index = find(bernstein_N_patches == 0,1,'first') + 1 ;
if isempty(problem_solved_index)
    problem_solved_index = default_max_iteration * 2 * dimension;
end
num_iter = floor((problem_solved_index) / (2*dimension)) ;

% get indices for all steps
all_indices = 1:problem_solved_index ;

% get the indices where elimination occurred
elim_indices = 2:2:problem_solved_index ;

% get the indices where subdivision occurred
subd_indices = all_indices ;
subd_indices(elim_indices) = [] ;

% get the maximum number of patches in each iteration
max_N_patches = max(reshape(bernstein_N_patches(:),2*dimension,[]),[],1) ;

% get x values for plotting the number of iterations
plot_subd_values = subd_indices ./ (2*dimension) ;
plot_elim_values = elim_indices ./ (2*dimension) ;

%% display problem info
disp(' ')
disp(['--- BENCHMARK P',num2str(problem_index),' STATS ---'])
disp(['Decision variable dimension: ',num2str(dimension)])
disp(['Iterations used: ',num2str(num_iter)])
disp(['Max number of patches: ',num2str(max(bernstein_N_patches))])
disp(['Max GPU memory used: ',num2str(max(bernstein_memory)),' ',mem_type,' (approx)'])
disp(' ')

%% plot
close all ; 
f1 = figure(1) ; clf ; hold on ;

% plot the max number of patches at every iteration
% plot(all_indices./(2*dimension),bernstein_N_patches(1:problem_solved_index),...
%     'b-','MarkerSize',12,'LineWidth',1);
plot(1:num_iter,max_N_patches(1:num_iter),...
    'b-','MarkerSize',12,'LineWidth',1);

% label left side
ylabel('max # of patches')

% add x ticks for iterations
grid on
xticks(1:2:num_iter+1)
xtickangle(45)

% get the yticks
yt_left = yticks ;

% plot memory usage on right
yyaxis right
% set(gca,'LineColor','k')
y_left_max = max(yt_left)*memory_per_item ;
h_mem = plot([0 1],[0,y_left_max],'LineStyle', 'none') ; % invisible plop
ylabel(['approximate memory used [',mem_type,']'])
set(gca,'Color','none')

% add title
title(['benchmark P',num2str(problem_index),': max # of patches / memory usage'])

% label plot
xlabel('iteration')
set(gca,'FontSize',14)

% set plot size
set(gcf,'Position',[500 400 600 400])

%% alternative plot
f2 = figure(2) ; clf ; hold on ;

% plot the number of patches at every subdivision step
h_subd = plot(plot_subd_values,bernstein_N_patches(subd_indices),'b.','MarkerSize',12) ;

% plot the number of patches at each elimination step
h_elim = plot(plot_elim_values,bernstein_N_patches(elim_indices),'r.','MarkerSize',12) ;

% label left side
ylabel('number of patches')

% add x ticks for iterations
grid on
xticks(1:2:num_iter)
xtickangle(45)

% get the yticks
yt_left = yticks ;

% plot memory usage on right
yyaxis right
% set(gca,'LineColor','k')
y_left_max = max(yt_left)*memory_per_item ;
h_mem = plot([0 1],[0,y_left_max],'LineStyle', 'none') ; % invisible plop
ylabel(['approximate memory used [',mem_type,']'])
set(gca,'Color','none')

% add legend
legend([h_subd h_elim],'subdivision','elimination')

% add title
title(['benchmark P',num2str(problem_index),' patches / memory usage'])

% label plot
xlabel('iteration')
set(gca,'FontSize',14)

% set plot size
set(gcf,'Position',[1000 800 600 400])

%% save figure
if save_pdf_flag
    % size the figure correctly
    set(f1,'Units','Inches');
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    % print
    print(f1,['P',num2str(problem_index),'_patches_vs_iterations.pdf'],'-dpdf','-r0')
end