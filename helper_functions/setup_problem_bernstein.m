function [bernstein_cost,bernstein_constraints,bernstein_cons_length] = setup_problem_bernstein(raw_cost,raw_constraints)
% setup problem input for bernstein mex function
%
% Author: Bohao Zhang
% Created: 6 Jan 2020

    bernstein_cost = raw_cost';

    bernstein_constraints = [];
    bernstein_cons_length = zeros(length(raw_constraints),1);
    for i = 1:length(raw_constraints)
        constraint = cell2mat(raw_constraints(i));
        bernstein_constraints = [bernstein_constraints;constraint];
        bernstein_cons_length(i) = size(constraint,1);
    end
    bernstein_constraints = bernstein_constraints';
end

