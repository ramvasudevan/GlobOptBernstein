function gradient = evaluate_function_gradient(problem_data,input)
% evaluate the gradient of the polynomial based on the data matrix of that function
%
% Author: Bohao Zhang
% Created: 6 Jan 2020

    numDimension = size(problem_data,2) - 1;
    numTerms = size(problem_data,1);
    gradient = zeros(1,numDimension);
    for dim = 1:numDimension
        value = 0;
        for i = 1:numTerms
            if problem_data(i,dim) > 0
                termValue = problem_data(i,numDimension + 1);
                for j = 1:numDimension
                    if j == dim
                        termValue = termValue * problem_data(i,j) * input(j) ^ (problem_data(i,j) - 1);
                    else
                        termValue = termValue * input(j) ^ problem_data(i,j);
                    end
                end
                value = value + termValue;
            end
        end
        gradient(dim) = value;
    end
end

