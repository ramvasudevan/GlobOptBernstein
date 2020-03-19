function [c, gc] = setup_cost_fmincon(cost_data,input)
% This function converts the cost function into fmincon format.
%
% Author: Bohao Zhang
% Created: 6 Jan 2020

    c = evaluate_function(cost_data,input);
    gc = evaluate_function_gradient(cost_data,input);
end

