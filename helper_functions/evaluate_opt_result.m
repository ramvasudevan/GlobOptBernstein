function [value,feasibility,violate_terms,difference] = evaluate_opt_result(cost,constraints,input)
% evaluate the value and the feasibility of the optimal solution found by
% the solver 
%
% Author: Bohao Zhang
% Created: 6 Jan 2020
violate_terms = [];
difference = [];
feasibility = 1;
for i = 1:length(constraints)
    cons_value = evaluate_function(cell2mat(constraints(i)),input);
    if cons_value > 0
        feasibility = 0;
        violate_terms = [violate_terms,i];
        if(cons_value > 1)
            difference = [difference,cons_value-1];
        else
            difference = [difference,cons_value];
        end
    end
end
value = evaluate_function(cost,input);
end
