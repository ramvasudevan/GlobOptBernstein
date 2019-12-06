function [c, gc] = setup_cost_fmincon(cost_data,input)
c = evaluate_function(cost_data,input);
gc = evaluate_function_gradient(cost_data,input);
end

