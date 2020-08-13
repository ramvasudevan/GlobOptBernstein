function [fmin,x] = DIRECT_increasing_constraints(numCons, problem_index, options)
% Increasing constraints problem setup for DIRECT
    load(['problem_matrix/more_',num2str(problem_index),'.mat']);
    
    bounds = repmat([0,1],size(more_problem.cost,2)-1,1);
    Problem.numconstraints = numCons;
    
    Problem.f = @(x)evaluate_function(more_problem.cost,x);
    
    for i = 1:numCons
        Problem.constraint(i).func = @(x)evaluate_function(more_problem.constraints{i},x);
        Problem.constraint(i).penalty = 1;
    end
    
    [fmin,x,~] = Direct(Problem,bounds,options);
end

