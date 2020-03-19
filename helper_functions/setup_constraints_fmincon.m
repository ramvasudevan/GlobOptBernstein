function [n, neq, gn, gneq] = setup_constraints_fmincon(constraints, equalities, input)
% This function converts constraints into fmincon format.
%
% Author: Bohao Zhang
% Created: 6 Jan 2020

    numDimension = length(input);
    numCons = length(constraints);
    numEqus = length(equalities);

    if numEqus > 0
        neq = zeros(numEqus,1);
        gneq = zeros(numDimension,numEqus);
        for i = 1:numEqus
            equ_data = cell2mat(equalities(i));
            neq(i) = evaluate_function(equ_data,input);
            gneq_ele = evaluate_function_gradient(equ_data,input);
            gneq(:,i) = gneq_ele';
        end
    else
        neq = [] ;
        gneq = [] ;
    end

    if numCons > 0
        % n = zeros(2*numCons,1);
        % gn = zeros(numDimension,2*numCons);
        % for i = 1:numCons
        %     constraint_data = cell2mat(constraints(i));
        %     n_ele = evaluate_function(constraint_data,input);
        %     n(2*i-1) = -n_ele;
        %     n(2*i) = n_ele - 1;
        %     gn_ele = evaluate_function_gradient(constraint_data,input);
        %     gn(:,2*i-1) = -gn_ele';
        %     gn(:,2*i) = gn_ele';
        % end
        n = zeros(numCons,1);
        gn = zeros(numDimension,numCons);
        for i = 1:numCons
            constraint_data = cell2mat(constraints(i));
            n(i) = evaluate_function(constraint_data,input);
            gn_ele = evaluate_function_gradient(constraint_data,input);
            gn(:,i) = gn_ele';
        end
    else
        n = [];
        gn = [];
    end
end

