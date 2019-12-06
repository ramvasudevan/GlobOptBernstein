function [F,G,I,J,d,k] = setup_problem_Lasserre(raw_cost,raw_constraints,d,k)
% setup problem for Lasserre
% k = ceil(max(max(raw_cost(:,1:(end-1)))) / 2); % the default k choice
% d = k; % the default d choice
F = raw_cost;
G = raw_constraints;
numDimension = size(raw_cost,2) - 1;

for i = 1:numDimension
    bound = zeros(1,numDimension+1);
    bound(numDimension+1) = 1;
    bound(i) = 1;
    G(end+1) = {bound};
end

I = {1:numDimension};
J = {1:length(G)};
end

