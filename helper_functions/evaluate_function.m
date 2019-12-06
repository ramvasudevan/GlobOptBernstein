function value = evaluate_function(problem_data,input)
% evaluate the function value according to the data matrix of that function
numDimension = size(problem_data,2) - 1;
numTerms = size(problem_data,1);
value = 0;
for i = 1:numTerms
    termValue = problem_data(i,numDimension + 1);
    for j = 1:numDimension
        termValue = termValue * (input(j) ^ problem_data(i,j));
    end
    value = value + termValue;
end
end

