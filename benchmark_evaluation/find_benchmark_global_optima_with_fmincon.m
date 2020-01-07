%% common parameters
% which problem(s) to solve
problem_indices = 3 ; % use 1:8 to run all problems (takes about 7 mins)

% solver params
opt = optimoptions('fmincon','Display','off',...
    'OptimalityTolerance',1e-10,...
    'MaxIterations',10^4,...
    'MaxFunctionEvaluations',10^5,...
    'ConstraintTolerance',1e-10) ;

% global min finding params
N_tries = 100 ;

% file i/o
save_file_flag = false ;
save_filename = 'PCBA_benchmark_optimal_solutions.mat' ;

%% automated from here
% set up to save info
optimal_value = cell(1,length(problem_indices)) ;
optimal_location = cell(1,length(problem_indices)) ;
time_spent = cell(1,length(problem_indices)) ;

total_time = tic ;
for problem_index = problem_indices
    [cost,nonlcon,lb,ub] = get_problem_data(problem_index) ;
    
    % set up best answer
    x = inf(length(lb),1) ;
    fval = inf ;
    
    time_in_loop = tic ;
    for idx = 1:N_tries
        % disp(['Attempt ',num2str(idx)])
        
        % new initial guess
        x0 = rand_range(lb-1,ub+1) ;
        
        % run fmincon
        [x_new,fval_new,exitflag] = fmincon(cost,x0,[],[],[],[],lb,ub,nonlcon,opt) ;
        
        % get best guess so far
        if (exitflag > 0) && (fval_new < fval)
            % disp('Answer improved!')
            x = x_new ;
            fval = fval_new ;
        end
    end
    time_in_loop = toc(time_in_loop) ;
    
    %% store output
    optimal_value{problem_index} = fval ;
    optimal_location{problem_index} = x ;
    time_spent{problem_index} = time_in_loop ;
    
end
total_time = toc(total_time) ;

%% display results
clc
for problem_index = problem_indices
    disp('--------')
    disp(['Problem ',num2str(problem_index)])
    disp(' ')
    
    disp('optimal value:')
    disp(['    ',num2str(optimal_value{problem_index},'%0.10f')])
    disp(' ')
    
    disp('optimal location:')
    for idx = 1:length(optimal_location{problem_index})
        disp(['    ',num2str(optimal_location{problem_index}(idx),'%0.10f')])
    end
end
disp(' ')
disp(['Total time spent: ',num2str(total_time/60),' min'])

%% save result
if save_file_flag
    save(save_filename,'problem_indices','optimal_value','optimal_location',...
        'time_spent','total_time') ;
end

%% problem costs and constraints
%% P1
function c = cost_1(x)
c = -x(1) - x(2) ;
end

function [n,neq] = nonlcon_1(x)
n = [-2*x(1)^4 + 8*x(1)^3 - 8*x(1)^2 + x(2) - 2 ;
    -4*x(1)^4 + 32*x(1)^3 - 88*x(1)^2 + 96*x(1) + x(2) - 36] ;
neq = [] ;
end

%% P2
function c = cost_2(x)
c = (x(1)-10)^3 + (x(2)-20)^3 ;
end

function [n,neq] = nonlcon_2(x)
n = [-(((x(1) - 5)^2) + ((x(2) - 5)^2) + 100) ;
    -(((x(1) - 5)^2) + ((x(2) - 5)^2) - 82.81)] ;
neq = [] ;
end

%% P3
function c = cost_3(x)
c = x(1) ;
end

function [n,neq] = nonlcon_3(x)
n = [x(1)^2 - x(2) ;
    x(2) - (x(1)^2)*(x(1) - 2) + (10^(-5))] ;
neq = [] ;
end

%% P4
function c = cost_4(x)
c = -2*x(1) + x(2) - x(3) ;
end

function [n,neq] = nonlcon_4(x)
A = [0 0 1 ; 0 -1 0 ; -2 1 -1] ;
b = [3;0;-4] ;
y = [1.5;-0.5;-5] ;
z = [0;-1;-6] ;

n = [-(x'*(A')*A*x - 2*(y')*A*x + norm(y)^2 - 0.25*(norm(b - z)^2)) ;
    x(1) + x(2) + x(3) - 4 ;
    3*x(2) + x(3) - 6] ;

neq = [] ;
end

%% P5
function c = cost_5(x)
c = x(3) ;
end

function [n,neq] = nonlcon_5(x)
neq = [] ;
n = [2*x(1)^2 + 4*x(1)*x(2) - 42*x(1) + 4*x(1)^3 - x(3) - 14 ;
    -2*x(1)^2 - 4*x(1)*x(2) + 42*x(1) - 4*x(1)^3 - x(3) + 14 ;
    2*x(1)^2 + 4*x(1)*x(2) - 26*x(2) + 4*x(2)^3 - x(3) - 22 ;
    -2*x(1)^2 - 4*x(1)*x(2) + 26*x(2) - 4*x(2)^3 - x(3) + 22 ] ;
end

%% P6
function c = cost_6(x)
c = 0.6224*x(3)*x(4) + 1.7781*x(2)*(x(3)^2) + 3.1661*(x(1)^2)*x(4) + 19.84*(x(1)^2)*x(3) ;
end

function [n,neq] = nonlcon_6(x)
neq = [] ;
n = [-x(1) + 0.0193*x(3) ;
    -x(2) + 0.00954*x(3) ;
    -pi*(x(3)^2)*x(4) - (4/3)*pi*(x(3)^3) + 750.1728 ;
    -240 + x(4)] ;
end

%% P7
function c = cost_7(x)
c = x(4) ;
end

function [n,neq] = nonlcon_7(x)
neq = (x(1)^4)*(x(2)^4) - (x(1)^4) - (x(2)^4)*x(3) ;
n = [1.4 - 0.25*x(4) - x(1) ;
    -1.4 - 0.25*x(4) + x(1) ;
    1.5 - 0.2*x(4) - x(2) ;
    -1.5 - 0.2*x(4) + x(2) ;
    0.8 - 0.2*x(4) - x(3) ;
    -0.8 - 0.2*x(4) + x(3) ] ;
end

%% P8
function c = cost_8(x)
c = 54.528*x(2)*x(4) + 27.264*x(1)*x(3) - 54.528*x(3)*x(4) ;
end

function [n,neq] = nonlcon_8(x)
neq = [] ;
I = 6*(x(1)^2)*x(2)*x(3) - 12*x(1)*x(2)*(x(3)^2) + 8*x(2)*(x(3)^3) + ...
    (x(1)^3)*x(4) - 6*(x(1)^2)*x(3)*x(4) + 12*x(1)*(x(3)^2)*x(4) - 8*(x(3)^3)*x(4) ;
n = [61.01627586 - I ;
    8*x(1) - I ;
    x(1)*x(2)*x(4) - x(2)*(x(4)^2) + (x(1)^2)*x(3) + x(3)*(x(4)^2) - 2*x(1)*x(3)*x(4) - 3.5*x(3)*I ;
    x(1) - 3*x(2) ;
    2*x(2) - x(1) ;
    x(3) - 1.5*x(4) ;
    0.5*x(4) - x(3)] ;
end

%% P2 - Bohao version
function c = cost_2_Bohao(x)
C = [0, 0, -7973;
    1, 0, 2349;
    2, 0, 68121;
    3, 0, 658500;
    0, 1, 120000;
    0, 2, -600000;
    0, 3, 1000000];
c = sum((x(1).^C(:,1).*x(2).^C(:,2).*C(:,3))) ;
end

function [n,neq] = nonlcon_2_Bohao(x)
neq = [] ;

G1 = [0, 0, 11;
    1, 0, -1392;
    2, 0, -7569;
    0, 1, 1000;
    0, 2, -10000];

G2 = [0, 0, -8.81;
    1, 0, 1218;
    2, 0, 7569;
    0, 1, -1000;
    0, 2, 10000];

n = [sum((x(1).^G1(:,1).*x(2).^G1(:,2).*G1(:,3))) ;
    sum((x(1).^G2(:,1).*x(2).^G2(:,2).*G2(:,3))) ] ;
end

%% helper functions
function [cost,nonlcon,lb,ub] = get_problem_data(problem_index)
switch problem_index
    case 1
        cost = @cost_1 ; nonlcon = @nonlcon_1 ;
        lb = [0 ; 0] ;
        ub = [3 ; 4] ;
    case 2
%         cost = @cost_2 ; nonlcon = @nonlcon_2 ;
%         lb = [13 ; 0] ;
%         ub = [100 ; 100] ;
        cost = @cost_2_Bohao ; nonlcon = @nonlcon_2_Bohao ;
        lb = [0 ; 0] ;
        ub = [1 ; 1] ;
    case 3
        cost = @cost_3 ; nonlcon = @nonlcon_3 ;
        lb = [-10 ; -10] ;
        ub = [10 ; 10] ;
    case 4
        cost = @cost_4 ; nonlcon = @nonlcon_4 ;
        lb = [0 ; 0 ; 0] ;
        ub = [2 ; 10 ; 3] ;
    case 5
        cost = @cost_5 ; nonlcon = @nonlcon_5 ;
        lb = [-5 ; -5 ; -5] ;
        ub = [5 ; 5 ; 5] ;
    case 6
        cost = @cost_6 ; nonlcon = @nonlcon_6 ;
        lb = [1 ; 0.625 ; 47.5 ; 90] ;
        ub = [1.375 ; 1 ; 52.5 ; 112] ;
    case 7
        cost = @cost_7 ; nonlcon = @nonlcon_7 ;
        lb = [0 ; 0 ; 0 ; 0] ;
        ub = [5 ; 5 ; 5 ; 5] ;
    case 8
        cost = @cost_8 ; nonlcon = @nonlcon_8 ;
        lb = [3 ; 2 ; 0.125 ; 0.25] ;
        ub = [20 ; 15 ; 0.75 ; 1.25] ;
    otherwise
        error('invalid problem index!')
end
end