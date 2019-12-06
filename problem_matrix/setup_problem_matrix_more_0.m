function [cost,constraints,feasible_point] = setup_problem_matrix_more_0(num)
% Cube
cost = [[   0,   0,  9801.01201]
        [   1,   0, -118800.044]
        [   2,   0,  597600.04]
        [   3,   0, -1598400]
        [   4,   0,  2400000]
        [   5,   0, -1920000]
        [   6,   0,  640000]
        [   0,   1,  396]
        [   0,   2,  4]
        [   1,   1, -2400]
        [   2,   1,  4800]
        [   3,   1, -3200]];

degree = [   0,     0;
             2,     0;
             0,     2;
             1,     1;
             1,     0;
             0,     1];
         
constraints = cell(num,1);

% This is the global optimum of Cube
% f(x) = 0
feasible_point = [0.55;
                  0.55];

for i = 1:num
    rand_num = 10 * rand(5,1) - 5;
    coef = [0;rand_num];
    con_mat = [degree,coef];
    diff = evaluate_function(con_mat,feasible_point);
    if diff >= 0
        con_mat(1,3) = con_mat(1,3) - diff - rand() / 10;
    end
    
    constraints(i) = {con_mat};
end

end



