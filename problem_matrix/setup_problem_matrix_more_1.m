function [cost,constraints,feasible_point] = setup_problem_matrix_more_1(num)
% Beale / 10000
cost = [[   0,   0,  10124.0741703125]
        [ 1.0,   0, -40486.6055]
        [ 2.0,   0,  40476.92]
        [   0, 1.0, -120946.16]
        [   0, 2.0,  602697.2]
        [   0, 3.0, -1603402]
        [   0, 4.0,  2401600]
        [   0, 5.0, -1920000]
        [   0, 6.0,  640000]
        [ 1.0, 1.0,  483725.12]
        [ 2.0, 1.0, -483665.6]
        [ 1.0, 2.0, -2410666.4]
        [ 2.0, 2.0,  2410544]
        [ 1.0, 3.0,  6413524]
        [ 2.0, 3.0, -6413440]
        [ 1.0, 4.0, -9606400]
        [ 2.0, 4.0,  9606400]
        [ 1.0, 5.0,  7680000]
        [ 2.0, 5.0, -7680000]
        [ 1.0, 6.0, -2560000]
        [ 2.0, 6.0,  2560000]];

degree = [   0,     0;
             2,     0;
             0,     2;
             1,     1;
             1,     0;
             0,     1];
         
constraints = cell(num,1);

% This is the global optimum of Beale / 10000
% f(x) = 0
feasible_point = [0.65;
                  0.525];

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