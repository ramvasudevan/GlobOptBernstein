% A script plotting a simple 1-D demonstration for the Bernstein algorithm
% with no constraints
%
% Author: Shreyas Kousik
% Created: 22 May 2021
clear ; clc ;
%% POP definition
% cost function coefficients
f_coeff = [470 -939  594 -124 4] ;

%% user parameters
% number of iterations to compute
N_iter = 8 ;

% plotting
plot_iterations = 1:N_iter  ;
plot_flag = true ;
plot_pause_duration = 0.25 ;

cost_color = [0 0.3 0.7] ;
cand_color = [0.8 0.8 1] ;
subo_color = [1 0.7 0.7] ;
optr_color = [0.2 0.5 0.1] ;

patch_linewidth = 1.5 ;
line_linewidth = 2 ;

f_bounds = [-6 6] ;

plot_position = [0 0 750 750] ;

plot_fontsize = 25 ;

% save figure output
save_png_flag = true ;

%% automated from here
% generate msspolys for f, g, and h
x = msspoly('x',1) ;

% make monomials for f and g, then generate f and g
x_f = monomials(x,0:(length(f_coeff) - 1)) ;

f = f_coeff(end:-1:1)*x_f ;

deg_f = deg(f) ;

%% list setup
% set up array to store the bernstein patches (each row is a patch)
B_f = zeros(1, deg_f + 3) ; % [center, width, coeffs]

% set up the width of the current box
dx = 1 ;

% current number of patches
N_patches = size(B_f,1) ;

% current solution estimate
solution_estimate = inf ;
optimizer_estimate = nan ;

% cell arrays for candidate optimal and definitely suboptimal patches
B_f_cand = cell(1,N_iter+1) ;
B_f_subo = cell(1,N_iter+1) ;

% this will store the row indices that are cut off in each iteration,
% whether they are infeasible or have a lower bound greater than the least
% upper bound
solution_estimate_all = cell(1,N_iter+1) ;
optimizer_estimate_all = cell(1,N_iter+1) ;

%% initialize
% set up the center and width of the first patch
B_f(1,1:2) = [0.5, 1] ; % center and width of first patch

% set up the bernstein coeffs (BCs) of the first patch
B_f(1,3:end) = Poly2Bernstein_1D(f) ;

% every patch is a candidate for the first iteration
B_f_cand{1} = B_f ;

% the initial solution estimate is infinity
solution_estimate_all{1} = solution_estimate ;
optimizer_estimate_all{1} = optimizer_estimate ;

%% run BA
for itr = 2:(N_iter+1)
    disp(['Running BA iteration ',num2str(itr-1)])
    
    %% begin current iteration
    % cut the patch width in half
    dx = dx / 2 ;
    
    % get previous iteration's candidate patches
    B_f = B_f_cand{itr-1} ;
    
    % get the new number of patches
    N_patches = size(B_f,1) ;
    
    % set up new patches to fill in
    B_f_new = zeros(2*N_patches, deg_f + 3) ;
    
    %% subdivide
    for idx = 1:size(B_f,1)
        % get center of current patch
        center_idx = B_f(idx,1) ;
        
        % get new BCs
        [f_coeff_L,f_coeff_R] = Bernstein_Dilation_1D(B_f(idx,3:end)) ;
        
        % fill in new patches
        B_f_new(2*idx-1,1) = center_idx - dx/2 ;
        B_f_new(2*idx-1,2) = dx ;
        B_f_new(2*idx-1,3:end) = f_coeff_L ;
        
        B_f_new(2*idx,1) = center_idx + dx/2 ;
        B_f_new(2*idx,2) = dx ;
        B_f_new(2*idx,3:end) = f_coeff_R ;
    end
    
    %% find bounds
    % extract the new patch bounds; we ignore the first two columns because
    % those contain the center and width of each patch, whereas the
    % remaining columns are the Bernstein coefficients
    B_f_bounds = [min(B_f_new(:,3:end),[],2), max(B_f_new(:,3:end),[],2)] ;
    
    %% cutoff test    
    % update solution estimate as smallest upper bound
    [solution_estimate,opt_row_idx] = min([solution_estimate ; B_f_bounds(:,2)]) ;
    solution_estimate_all{itr} = solution_estimate ;
    
    % update optimizer as center of patch with smallest upper bound, or
    % leave unchanged
    if opt_row_idx > 1
        optimizer_estimate = B_f_new(opt_row_idx-1,1) ;
    end
    optimizer_estimate_all{itr} = optimizer_estimate ;
    
    % get suboptimal patches
    subo_log = (B_f_bounds(:,1) > solution_estimate) ; 
    
    %% eliminate
    % save candidate patches
    B_f_cand{itr} = B_f_new(~subo_log,:) ;
    
    % save suboptimal patches
    B_f_subo{itr} = B_f_new(subo_log,:) ;
end

%% plotting
if plot_flag
    % prep for plotting
    x_vals = 0:0.001:1 ;
    f_vals = polyval(f_coeff,x_vals) ;
    
    ccolor = cand_color ;
    scolor = subo_color ;
    lw = patch_linewidth ;
    
    close all
    
    for itr = plot_iterations
        disp(['Plotting BA iteration ',num2str(itr)])
        
        %% plot
        fh = figure(1) ; clf ; hold on ; axis([0 1 f_bounds])
        if itr == 1
            set(fh,'Position',plot_position) ;
        end

        % plot patches
        C = B_f_cand{itr} ;
        S = B_f_subo{itr} ;
        plot_patches(C,S,ccolor,scolor,lw)
        
        % plot cost
        plot(x_vals,f_vals,'Color',cost_color,'LineWidth',line_linewidth) ;
        
        % plot solution estimate
        plot([0 1],solution_estimate_all{itr}.*ones(1,2),'--',...
            'Color',optr_color,'LineWidth',line_linewidth) ;
        
        % plot optimizer estimate
        plot(optimizer_estimate_all{itr}.*ones(1,2),f_bounds,'--',...
            'Color',optr_color,'LineWidth',line_linewidth) ;
        
        % labels and stuff
        xlabel('decision variable')
        ylabel('cost')
        set(gca,'FontSize',plot_fontsize,'LineWidth',1)
        
        % plot and pause for dramatic effect
        drawnow
        pause(plot_pause_duration)
        
        %% save output
        if save_png_flag
            save_png_filename = ['ba_no_cons_iter_',num2str(itr),'.png'] ;
            save_figure_to_png(fh,save_png_filename) ;
        end
    end
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

function plot_patch(patch_info,color,varargin)
    c = patch_info(1) ; % center
    w = patch_info(2) ; % width
    l = min(patch_info(3:end)) ; % lower bound
    u = max(patch_info(3:end)) ; % upper bound

    % make vertices
    x = [c-w/2, c+w/2, c+w/2, c-w/2] ;
    y = [l,l,u,u] ;

    % plot patch
    patch(x,y,color,varargin{:}) ;

    % % plot an x in a patch if it is to be eliminated
    % if elim_flag
    %     x = c ;
    %     y = (l+u)/2 ;
    %     plot(x,y,'kx','LineWidth',2,'MarkerSize',10)
    % end
end

function plot_patches(C,S,ccolor,scolor,lw)
    % plot candidate patches
    for idx = 1:size(C,1)
        plot_patch(C(idx,:),ccolor,'LineWidth',lw)
    end

    % plot suboptimal patches
    for idx = 1:size(S,1)
        plot_patch(S(idx,:),scolor,'LineWidth',lw)
    end
end