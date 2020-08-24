%% user parameters
% iteration to plot
iter = 11 ; % pick in 1:49

% colors
cost_line_color = [0 0.3 1] ;
cons_line_color = [0.5 0 0.5] ;

cost_fill_color = [0.3 0.6 1] ;
cons_fill_color = [1 0.7 1] ;

feas_area_color = [1 1 1] ;

% save output
save_pdfs_flag = false ;

%% automated from here
% load the data
load('demo_2_PCBA_POP_data.mat')
load('segway_FRSes.mat')

%% get the cost and constraint data
cdata = cdata_cell{iter} ;
cons_poly = cons_poly_cell{iter} ;

% set up msspoly
q = msspoly('q',2) ;

% reconstitute the cost and constraints from the data
p = recomp(q,cdata(:,1:2),cdata(:,3)') ;
g = recomp(q,cons_poly.pows,cons_poly.coef) ;

% make grid for plotting
q_vec = linspace(0,1) ;
[Q1,Q2] = meshgrid(q_vec,q_vec) ;
Q = [Q1(:) Q2(:)]' ;

% evaluate cost and constraints on grid
P = reshape(full(msubs(p,q,Q)),100,100) ;
G = cell(1,length(g)) ;

tic
for idx = 1:length(g)
    G{idx} = reshape(full(msubs(g(idx),q,Q)),100,100) ;
end
toc

%% prep to plot agent in world
% get segway info
z0 = z0_cell{iter} ;

% make circle for segway plot
footprint = 0.38 ;
t_body = linspace(0,2*pi) ;
body_vertices = [footprint.*cos(t_body) + z0(1) ;
                 footprint.*sin(t_body) + z0(2) ] ;
                 
% make arrow for plot
t_arrow = [0, 2*pi/3, 4*pi/3] ;
x_arrow = 0.4*footprint*cos(t_arrow) ;
y_arrow = 0.3*footprint*sin(t_arrow) ;
R = rotation_matrix_2D(z0(3)) ;
arrow_vertices = R*[x_arrow ; y_arrow] + repmat(z0(1:2),1,3) ;

% set up patch data for segway
body.faces = [1:(size(body_vertices,2)-1), 1] ;
body.vertices = body_vertices' ;
arrow.faces = [1 2 3 1] ;
arrow.vertices = arrow_vertices' ;

% get obstacle data
O = O_static_cell{iter} ;

%%% HACK HERE FOR REMOVING ONE STRAY POINT
O_log = (O(1,:) > 0) & (O(1,:) < 1) & (O(2,:) > 0) ;
O(:,O_log) = [] ;
%%% END HACK FOR REMOVING ONE STRAY POINT

% get waypoint
wp = wp_cell{iter} ;

% move waypoint to global coordinates
wp = local_to_world(z0,wp) ;

%% create trajectory to plot in world
% find which FRS was used
Dx = Dx_cell{iter} ;
Dy = Dy_cell{iter} ;
x0 = x0_cell{iter} ;
y0 = y0_cell{iter} ;

FRS_idx = nan ;

for idx = 1:length(FRSes)
    Dx_idx = FRSes{idx}.Dx_static ;
    Dy_idx = FRSes{idx}.Dy_static ;
    x0_idx = FRSes{idx}.x0_static ;
    y0_idx = FRSes{idx}.y0_static ;
    
    if Dx_idx == Dx && Dy_idx == Dy && x0_idx == x0 && y0_idx == y0
        FRS = FRSes{idx} ;
        break
    end
end

%%
% get optimal solution
q_opt = k_opt_pcba_cell{iter} ;

% scale q opt to [-1,1]
q_scaled = 2*(q_opt - 0.5) ;

% convert q opt into velocity and yaw rate commands
v_min = FRS.min_velocity ;
v_max = FRS.max_velocity ;
w_min = FRS.min_yawrate ;
w_max = FRS.max_yawrate ;

v_des = ((v_max - v_min)/2)*q_scaled(2) + ((v_max + v_min)/2) ;
w_des = ((w_max - w_min)/2)*q_scaled(2) + ((w_max + w_min)/2) ;

% create trajectory for plotting
u = [w_des ; v_des] ;
[~,z_des] = ode45(@(t,z) traj_dyn(t,z,u),[0 FRS.T_static],z0(1:3)) ;
z_des = z_des(:,1:2)' ;

%% plot cost and constraints
f_q = figure(1) ; clf ; hold on ; axis equal

% background patch
patch('Faces',[1 2 3 4 1],'Vertices',[0 1 1 0 ; 0 0 1 1]','FaceColor',feas_area_color)

% constraint contours
colormap([cons_fill_color ; cons_fill_color])
for idx = 1:length(G)
    contourf(q_vec,q_vec,G{idx},[0 0],'LineWidth',1,'LineColor',cons_line_color) ;
end

% cost contour
c_p = contour(Q1,Q2,2.*P - 9,-10:20,'LineColor',cost_line_color,'LineWidth',1.5) ;
clabel(c_p,'FontSize',14,'Color',cost_line_color,'FontWeight','bold')

% optimal solution
plot(q_opt(1),q_opt(2),'p','MarkerFaceColor',cost_line_color,'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',16)

% labels
xlabel('q_1')
ylabel('q_2')

set(gca,'FontSize',14)

%% plot world
f_z = figure(2) ; clf ; hold on ; axis equal ;

plot_path(O,'.','Color',cons_line_color,'MarkerSize',10)
    plot(wp(1),wp(2),'p','MarkerFaceColor',cost_line_color,'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',16)
patch(body,'FaceColor',[0 0 1],'EdgeColor',[0 0 0],'FaceAlpha',0.1)
patch(arrow,'FaceColor',[0 0 1],'EdgeColor',[0 0 0],'FaceAlpha',0.5)
plot_path(z_des,'--','Color',cost_line_color,'LineWidth',1.5)

%%% custom axes for the plotted iteration
axis([-1.5 2 -3.5 0.6])
xlabel('x [m]')
ylabel('y [m]')

set(gca,'FontSize',14)

%% save output
if save_pdfs_flag
    save_figure_to_pdf(f_q,'Q_space_cost_and_constraints.pdf')
    save_figure_to_pdf(f_z,'Z_space_obs_robot_and_wp.pdf')
end

%% helper function
function P_out = local_to_world(robot_pose, P_local)
    % extract position and heading from input
    x = robot_pose(1,1) ;
    y = robot_pose(2,1) ;
    h = robot_pose(3,1) ;
    
    % prep
    P_out = P_local ;
    [N_rows,N_cols] = size(P_local) ;
    
    % rotate points to world frame direction
    R = [cos(h), -sin(h) ;
         sin(h),  cos(h) ] ;
    P_out(1:2,:) = R*P_out(1:2,:) ;
    
    if N_rows > 2
        P_out(3,:) = P_out(3,:) + h ;
    end
    
    % shift points to world frame location
    P_out(1:2,:) = P_out(1:2,:) + repmat([x;y],1,N_cols) ;
end

function dz = traj_dyn(~,z,u)
    dz = [u(2)*cos(z(3)) ;
        u(2)*sin(z(3)) ;
        u(1)] ;
end