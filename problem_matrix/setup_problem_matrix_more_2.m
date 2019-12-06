function [cost,constraints,feasible_point] = setup_problem_matrix_more_2(num)
% Bukin02
cost = [[   0,   0, 775.25]
        [ 1.0,   0, 299   ]
        [ 2.0,   0, -99   ]
        [   0, 1.0, -3600 ]
        [   0, 2.0, 3600  ]];

degree = [   0,     0;
             2,     0;
             0,     2;
             1,     1;
             1,     0;
             0,     1];
         
constraints = cell(num,1);

% This is the global optimum of Bukin02
% f(x) = -124.75
feasible_point = [0;
                  0.5];

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