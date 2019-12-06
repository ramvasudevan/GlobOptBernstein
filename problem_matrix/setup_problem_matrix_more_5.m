function [cost,constraints,feasible_point] = setup_problem_matrix_more_5(num)
% DixonPrice 3D / 1000
cost = [[   0,   0,   0,                           220.621]
        [ 1.0,   0,   0,                            -17.24]
        [ 2.0,   0,   0,                               1.2]
        [   0, 1.0,   0,                            -697.2]
        [   0, 2.0,   0,                            1953.2]
        [   0, 3.0,   0,                           -2560.0]
        [   0, 4.0,   0,                            1280.0]
        [ 1.0, 1.0,   0,                              64.0]
        [ 1.0, 2.0,   0,                             -64.0]
        [   0,   0, 1.0,                           -1008.0]
        [   0,   0, 2.0,                            2928.0]
        [   0,   0, 3.0,                           -3840.0]
        [   0,   0, 4.0,                            1920.0]
        [   0, 1.0, 1.0,                              96.0]
        [   0, 1.0, 2.0,                             -96.0]];

degree = [  0, 0, 0;
            2, 0, 0;
            0, 2, 0;
            0, 0, 2;
            1, 1, 0;
            0, 1, 1;
            1, 0, 1;
            1, 0, 0;
            0, 1, 0;
            0, 0, 1];
         
constraints = cell(num,1);

% This is the global optimum of DixonPrice 3D / 1000
% f(x) = 0
feasible_point = ([2^(-(2^1-2)/(2^1));
                   2^(-(2^2-2)/(2^2));
                   2^(-(2^3-2)/(2^3))] + 10) ./ 20;

for i = 1:num
    rand_num = 10 * rand(9,1) - 5;
    coef = [0;rand_num];
    con_mat = [degree,coef];
    diff = evaluate_function(con_mat,feasible_point);
    if diff >= 0
        con_mat(1,4) = con_mat(1,4) - diff - rand() / 10;
    end
    
    constraints(i) = {con_mat};
end

end

