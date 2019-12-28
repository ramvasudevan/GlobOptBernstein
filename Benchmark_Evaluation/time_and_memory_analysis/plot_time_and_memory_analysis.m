%% description
% This script plots the number of Bernstein patches (i.e., the length of
% the list) versus the number of iterations of PCBA for each of the 8
% benchmark problems.
%
% Author: Shreyas Kousik
% Created: 28 Dec 2019
% Updated: -
%
%% user parameters
% which problem to plot
problem_index = 4 ;

% whether or not to save the output
save_pdf_flag = false ;

%% automated from here
% load data
load(['P',num2str(problem_index),'_time_and_memory_info.mat'])

%% process data
% get the problem dimension
d = size(bernstein_cost,1) - 1 ;

% get the number of iterations (each entry in bernstein_memory is the
% number of patches either from subdivision or elimination, and subdivision
% happens once along each dimension in each iteration, so the denominator
% is 2*d; the first "0" entry in bernstein_memory marks the iteration
% where the problem was solved to the specified tolerances)
problem_solved_index = find(bernstein_memory == 0,1,'first') + 1 ;
num_iter = (problem_solved_index) / (2*d) ;

% get indices for all steps
all_indices = 1:problem_solved_index ;

% get the indices where elimination occurred
elim_indices = 2:2:problem_solved_index ;

% get the indices where subdivision occurred
subd_indices = all_indices ;
subd_indices(elim_indices) = [] ;

% get x values for plotting the number of iterations
plot_subd_values = subd_indices ./ (2*d) ;
plot_elim_values = elim_indices ./ (2*d) ;

%% plot
close all ; 
f1 = figure(1) ; clf ; hold on ;

% plot the number of patches at every subdivision step
h_subd = plot(plot_subd_values,bernstein_memory(subd_indices),'b.','MarkerSize',12) ;

% plot the number of patches at each elimination step
h_elim = plot(plot_elim_values,bernstein_memory(elim_indices),'r.','MarkerSize',12) ;

% add x ticks for every iteration
xticks(1:num_iter)
xtickangle(45)

% add legend
legend([h_subd h_elim],'Subdivision','Elimination','Location','NorthWest')

% label plot
xlabel('Iteration')
ylabel('Number of Patches')
grid on
set(gca,'FontSize',14)

% set plot size
set(gcf,'Position',[1000 800 873 291])

%% save figure
if save_pdf_flag
    % size the figure correctly
    set(f1,'Units','Inches');
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    % print
    print(f1,['P',num2str(problem_index),'_patches_vs_iterations.pdf'],'-dpdf','-r0')
end

