function new_constraints = scale_for_lss(constraints)
% scale the constraint to [0,1]
numCons = length(constraints);
numDimension = size(cell2mat(constraints(1)),2) - 1;
new_constraints = cell(numCons,1);
for conID = 1:numCons 
    constraint_mat = cell2mat(constraints(conID));
    numTerms = size(constraint_mat,1);
    degrees = max(constraint_mat(:,1:numDimension)) + 1;
    con_unitLength = prod(degrees);
    con_BC = zeros(con_unitLength,1);
    for i = 1:numTerms
        for j = 1:con_unitLength
            index = j-1;
            BCterm = 1;
            for k = 1:numDimension
                currentDegree = mod(index, degrees(k));
                if currentDegree < constraint_mat(i,k)
                    BCterm = 0;
                else
                    BCterm = BCterm * nchoosek(currentDegree,constraint_mat(i,k)) / nchoosek(degrees(k)-1,constraint_mat(i,k));
                end
                index = (index - currentDegree) / degrees(k);
            end
            con_BC(j) = con_BC(j) + constraint_mat(i,numDimension+1) * BCterm;
        end
    end
    
    upp = max(con_BC);
    downn = min(con_BC);
    
%     x = msspoly('x',numDimension);
%     original_poly = recomp(x,constraint_mat(:,1:numDimension),constraint_mat(:,numDimension+1)');
%     if downn < -1
%         new_poly = original_poly / abs(downn) + 1;
%     else
%         new_poly = original_poly + 1;
%     end
%     [~,cpows,ccoef] = decomp(new_poly,x);
%     new_cons_degree = full(cpows);
%     new_cons_coef = full(ccoef);
%     new_constraints(conID) = {[new_cons_degree,new_cons_coef']};
    
    if downn < -1
        constraint_mat(:,numDimension+1) = constraint_mat(:,numDimension+1) / abs(downn);
        constraint_mat(1,numDimension+1) = constraint_mat(1,numDimension+1) + 1;
    else
        constraint_mat(1,numDimension+1) = constraint_mat(1,numDimension+1) + 1;
    end
    new_constraints(conID) = {constraint_mat};
end

end

