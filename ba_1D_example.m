% A script plotting a simple 1-D demonstration for the Bernstein algorithm
%
% Author: Shreyas Kousik
% Created: 6 Jan 2020

%% POP definition
% cost function coefficients
f_coeff = [470 -939  594 -124 4] ;

% constraint function coefficients
g_coeff = [-76, 118, -45, 3] ;

% equality constraint coefficients
h_coeff = 2 * 1.0e+05 * [0.1170
   -0.1961
   -0.5965
    2.1800
   -2.8760
    2.0204
   -0.8119
    0.1831
   -0.0209
    0.0009
   -0.0000]' ;

%% user parameters
% number of iterations to compute
N_iter = 10 ;

% equality constraint tolerance
eeq = 1 ;

% plotting
plot_iterations = 4  ;
plot_flag = true ;
plot_pause_duration = 0.25 ;

cost_color = [0 0.3 0.7] ;
cons_color = [0.5 0 0.5] ;

feas_color = [0.7 1 0.7] ;
infs_color = [1 0.7 1] ;
undc_color = [0.8 0.8 0.8] ;
subo_color = [0.7 1 1] ;

patch_linewidth = 1 ;
line_linewidth = 2 ;

f_bounds = [-6 6] ;
g_bounds = [-4 4] ;
h_bounds = [-3 3] ;

plot_position = [0 0 750 750] ;

% save figure output
save_pdf_flag = true ;
save_pdf_filename = ['pcba_demo_iteration_',num2str(plot_iterations(end)),'.pdf'] ;

%% automated from here
% generate msspolys for f, g, and h
x = msspoly('x',1) ;

% make monomials for f and g, then generate f and g
x_f = monomials(x,0:(length(f_coeff) - 1)) ;
x_g = monomials(x,0:(length(g_coeff) - 1)) ;
x_h = monomials(x,0:(length(h_coeff) - 1)) ;

f = f_coeff(end:-1:1)*x_f ;
g = g_coeff(end:-1:1)*x_g ;
h = h_coeff(end:-1:1)*x_h ;

deg_f = deg(f) ;
deg_g = deg(g) ;
deg_h = deg(h) ;

%% list setup
% set up array to store the bernstein patches (each row is a patch)
B_f = zeros(1, deg_f + 3) ; % [center, width, coeffs]
B_g = zeros(1, deg_g + 3) ;
B_h = zeros(1, deg_h + 3) ;

% set up the width of the current box
dx = 1 ;

% current number of patches
N_patches = size(B_f,1) ;

% current solution estimate
solution_estimate = inf ;

% cell arrays for feasible, infeasible, undecided, and suboptimal patches
B_f_feas = cell(1,N_iter+1) ;
B_g_feas = cell(1,N_iter+1) ;
B_h_feas = cell(1,N_iter+1) ;

B_f_infs = cell(1,N_iter+1) ;
B_g_infs = cell(1,N_iter+1) ;
B_h_infs = cell(1,N_iter+1) ;

B_f_undc = cell(1,N_iter+1) ;
B_g_undc = cell(1,N_iter+1) ;
B_h_undc = cell(1,N_iter+1) ;

B_f_subo = cell(1,N_iter+1) ;
B_g_subo = cell(1,N_iter+1) ;
B_h_subo = cell(1,N_iter+1) ;

% this will store the row indices that are cut off in each iteration,
% whether they are infeasible or have a lower bound greater than the least
% upper bound
solution_estimate_all = cell(1,N_iter+1) ;

%% initialize
% set up the center and width of the first patch
B_f(1,1:2) = [0.5, 1] ; % center and width of first patch
B_g(1,1:2) = [0.5, 1] ;
B_h(1,1:2) = [0.5, 1] ;

% set up the bernstein coeffs (BCs) of the first patch
B_f(1,3:end) = Poly2Bernstein_1D(f) ;
B_g(1,3:end) = Poly2Bernstein_1D(g) ;
B_h(1,3:end) = Poly2Bernstein_1D(h) ;

% every patch is undecided for the first iteration
B_f_undc{1} = B_f ;
B_g_undc{1} = B_g ;
B_h_undc{1} = B_h ;

% the initial solution estimate is infinity
solution_estimate_all{1} = solution_estimate ;

%% run BA
for itr = 2:(N_iter+1)
    %% begin current iteration
    % cut the patch width in half
    dx = dx / 2 ;
    
    % get previous iteration's patches
    B_f = [B_f_undc{itr-1}; B_f_feas{itr-1}] ;
    B_g = [B_g_undc{itr-1}; B_g_feas{itr-1}] ;
    B_h = [B_h_undc{itr-1}; B_h_feas{itr-1}] ;
    
    % get the new number of patches
    N_patches = size(B_f,1) ;
    
    % set up new patches to fill in
    B_f_new = zeros(2*N_patches, deg_f + 3) ;
    B_g_new = zeros(2*N_patches, deg_g + 3) ;
    B_h_new = zeros(2*N_patches, deg_h + 3) ;
    
    %% subdivide
    for idx = 1:size(B_f,1)
        % get center of current patch
        center_idx = B_f(idx,1) ;
        
        % get new BCs
        [f_coeff_L,f_coeff_R] = Bernstein_Dilation_1D(B_f(idx,3:end)) ;
        [g_coeff_L,g_coeff_R] = Bernstein_Dilation_1D(B_g(idx,3:end)) ;
        [h_coeff_L,h_coeff_R] = Bernstein_Dilation_1D(B_h(idx,3:end)) ;
        
        % fill in new patches
        B_f_new(2*idx-1,1) = center_idx - dx/2 ;
        B_f_new(2*idx-1,2) = dx ;
        B_f_new(2*idx-1,3:end) = f_coeff_L ;
        
        B_f_new(2*idx,1) = center_idx + dx/2 ;
        B_f_new(2*idx,2) = dx ;
        B_f_new(2*idx,3:end) = f_coeff_R ;
        
        B_g_new(2*idx-1,1) = center_idx - dx/2 ;
        B_g_new(2*idx-1,2) = dx ;
        B_g_new(2*idx-1,3:end) = g_coeff_L ;
        
        B_g_new(2*idx,1) = center_idx + dx/2 ;
        B_g_new(2*idx,2) = dx ;
        B_g_new(2*idx,3:end) = g_coeff_R ;
        
        B_h_new(2*idx-1,1) = center_idx - dx/2 ;
        B_h_new(2*idx-1,2) = dx ;
        B_h_new(2*idx-1,3:end) = h_coeff_L ;
        
        B_h_new(2*idx,1) = center_idx + dx/2 ;
        B_h_new(2*idx,2) = dx ;
        B_h_new(2*idx,3:end) = h_coeff_R ;
    end
    
    %% find bounds
    B_f_bounds = [min(B_f_new(:,3:end),[],2), max(B_f_new(:,3:end),[],2)] ;
    B_g_bounds = [min(B_g_new(:,3:end),[],2), max(B_g_new(:,3:end),[],2)] ;
    B_h_bounds = [min(B_h_new(:,3:end),[],2), max(B_h_new(:,3:end),[],2)] ;
    
    %% cutoff test
    % get feasible patches
    eq_cons_log = (B_h_bounds(:,1) >= (-eeq)) & (B_h_bounds(:,2) <= eeq) & ...
        (B_h_bounds(:,1) <= 0) & (B_h_bounds(:,2) >= 0) ;
    ineq_cons_log = B_g_bounds(:,2) <= 0 ;
    feas_log = eq_cons_log & ineq_cons_log ;
    
    % update solution estimate as smallest upper bound of feasible patches
    solution_estimate = min([solution_estimate ; B_f_bounds(feas_log,2)]) ;
    solution_estimate_all{itr} = solution_estimate ;
    
    % get infeasible patches
    infs_log = (B_h_bounds(:,1) > eeq) | (B_h_bounds(:,2) < -eeq) | ...
        (B_h_bounds(:,1) > 0) | (B_h_bounds(:,2) < 0) | (B_g_bounds(:,1) > 0) ;
    
    % get undecided patches
    undc_log = (~infs_log) & (~feas_log) ;
    
    % get suboptimal patches
    subo_log = (B_f_bounds(:,1) > solution_estimate) & (undc_log | feas_log) ;
    
    %% eliminate
    % save undecided and feasible patches
    B_f_undc{itr} = B_f_new(undc_log,:) ;
    B_g_undc{itr} = B_g_new(undc_log,:) ;
    B_h_undc{itr} = B_h_new(undc_log,:) ;
    
    B_f_feas{itr} = B_f_new(feas_log,:) ;
    B_g_feas{itr} = B_g_new(feas_log,:) ;
    B_h_feas{itr} = B_h_new(feas_log,:) ;
    
    % eliminate infeasible and suboptimal patches
    B_f_infs{itr} = B_f_new(infs_log,:) ;
    B_g_infs{itr} = B_g_new(infs_log,:) ;
    B_h_infs{itr} = B_h_new(infs_log,:) ;
    
    B_f_subo{itr} = B_f_new(subo_log,:) ;
    B_g_subo{itr} = B_g_new(subo_log,:) ;
    B_h_subo{itr} = B_h_new(subo_log,:) ;
end

%% plotting
if plot_flag
    % prep for plotting
    x_vals = 0:0.001:1 ;
    f_vals = polyval(f_coeff,x_vals) ;
    g_vals = polyval(g_coeff,x_vals) ;
    h_vals = polyval(h_coeff,x_vals) ;
    
    fcolor = feas_color ;
    icolor = infs_color ;
    ucolor = undc_color ;
    scolor = subo_color ;
    lw = patch_linewidth ;
    
    close all
    
    for itr = plot_iterations
        fh = figure(1) ; clf ; hold on ;
        set(fh,'Position',plot_position) ;
        
        %% plot cost
        subplot(3,1,1) ; hold on ; axis([0 1 f_bounds])

        % plot patches
        F = B_f_feas{itr} ;
        I = B_f_infs{itr} ;
        U = B_f_undc{itr} ;
        S = B_f_subo{itr} ;
        plot_patches(F,I,U,S,fcolor,icolor,ucolor,scolor,lw)
        
        % plot cost
        plot(x_vals,f_vals,'Color',cost_color,'LineWidth',line_linewidth) ;
        
        % plot solution estimate
        plot([0 1],solution_estimate_all{itr}.*ones(1,2),'--',...
            'Color',cost_color,'LineWidth',line_linewidth) ;
        
        % labels and stuff
        %xlabel('decision variable')
        ylabel('cost')
        %ylabel('$p(x)$','Interpreter','latex')
        set(gca,'FontSize',17,'LineWidth',1.5)
        
        %% plot inequality constraint
        subplot(3,1,2) ; hold on ; axis([0 1 g_bounds])

        % plot patches
        F = B_g_feas{itr} ;
        I = B_g_infs{itr} ;
        U = B_g_undc{itr} ;
        S = B_g_subo{itr} ;
        plot_patches(F,I,U,S,fcolor,icolor,ucolor,scolor,lw)
        
        % plot inequality constraint
        plot(x_vals,g_vals,'Color',cons_color,'LineWidth',line_linewidth) ;
        
        % plot 0 line
        plot([0 1],[0 0],'--','Color',cons_color,'LineWidth',line_linewidth) ;
        
        % labels and stuff
        %xlabel('decision variable')
        ylabel('ineq. cons.')
        %ylabel('$g(x)$','Interpreter','latex')
        set(gca,'FontSize',17,'LineWidth',1.5)
        
        %% plot equality constraint
        subplot(3,1,3) ; hold on ; axis([0 1 h_bounds])

        % plot patches
        F = B_h_feas{itr} ;
        I = B_h_infs{itr} ;
        U = B_h_undc{itr} ;
        S = B_h_subo{itr} ;
        plot_patches(F,I,U,S,fcolor,icolor,ucolor,scolor,lw)
        
        % plot equality constraint
        plot(x_vals,h_vals,'Color',cons_color,'LineWidth',line_linewidth) ;
        
        % plot 0 line
        plot([0 1],[0 0],'-','Color','k','LineWidth',line_linewidth/1.75) ;
        
        % plot equality constraint tolerances
        plot([0 1],[eeq eeq],'--','Color',cons_color,'LineWidth',line_linewidth) ;
        plot([0 1],[-eeq -eeq],'--','Color',cons_color,'LineWidth',line_linewidth) ;
        
        % labels and stuff
        xlabel('decision variable')
        ylabel('eq. cons.')
        %xlabel('decision variable $x$','Interpreter','latex')
        %ylabel('$h(x)$','Interpreter','latex')
        set(gca,'FontSize',17,'LineWidth',1.5)
        
        %% plot and pause for dramatic effect
        drawnow
        pause(plot_pause_duration)
    end
end

%% save output
if save_pdf_flag
    save_figure_to_pdf(fh,save_pdf_filename) ;
end

%% helper functions
function [ coeff ] = Poly2Bernstein_1D( p, n )
[ ~, pow, M ] = decomp( p );
if nargin < 2
    n = deg( p );
end
if deg( p ) == 0
    coeff = ones( n+1, 1 );
    return;
end
a = zeros( n + 1, 1 );
a( pow + 1 ) = M;

Mat = zeros( n+1, n+1 );

for i = 0 : n
    for j = 0 : i
        Mat( i+1, j+1 ) = nchoosek( i, j ) / nchoosek( n, j );
    end
end

coeff = Mat * a;
end

function [ coeffA, coeffB ] = Bernstein_Dilation_1D( coeff )

n = length( coeff ) - 1;
coeffA = coeff;
coeffB = zeros( 1, n+1 );
coeffB( end ) = coeffA( end );

% This is from Nataraj and Arounassalame 2007 Equation (5), with
% \lambda = 0.5
for k = 1 : n
    tmpA = coeffA;
    for i = k : n
        coeffA( i + 1 ) = 0.5 * ( tmpA(i) + tmpA(i+1) );
    end
    coeffB( n - k + 1 ) = coeffA( end );
end
% coeffB = flipud( coeffA );
end

function plot_patch(patch_info,color,elim_flag,varargin)
c = patch_info(1) ; % center
w = patch_info(2) ; % width
l = min(patch_info(3:end)) ; % lower bound
u = max(patch_info(3:end)) ; % upper bound

% make vertices
x = [c-w/2, c+w/2, c+w/2, c-w/2] ;
y = [l,l,u,u] ;

% plot patch
patch(x,y,color,varargin{:}) ;

% plot an x in a patch if it is to be eliminated
if elim_flag
    x = c ;
    y = (l+u)/2 ;
    plot(x,y,'kx','LineWidth',2,'MarkerSize',10)
end
end

function plot_patches(F,I,U,S,fcolor,icolor,ucolor,scolor,lw)
% plot feasible patches
for idx = 1:size(F,1)
    plot_patch(F(idx,:),fcolor,false,'LineWidth',lw)
end

% plot infeasible patches
for idx = 1:size(I,1)
    plot_patch(I(idx,:),icolor,true,'LineWidth',lw)
end

% plot undecided patches
for idx = 1:size(U,1)
    plot_patch(U(idx,:),ucolor,false,'LineWidth',lw)
end

% plot suboptimal patches
for idx = 1:size(S,1)
    plot_patch(S(idx,:),scolor,true,'LineWidth',lw)
end
end