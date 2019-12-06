function [cost,constraints,feasible_point] = setup_problem_matrix_more_7(num)
% El-Attar-Vidyasagar-Dutta
cost = [[   0,   0,                       980297.66355]
        [ 1.0,   0,                         78412.8372]
        [ 2.0,   0,                          -75208.84]
        [ 3.0,   0,                            -6400.0]
        [ 4.0,   0,                             3200.0]
        [   0, 1.0,                      -11880799.484]
        [   0, 2.0,                        59762415.48]
        [   0, 3.0,                       -159843216.0]
        [   0, 4.0,                        240001600.0]
        [   0, 5.0,                       -192000000.0]
        [   0, 6.0,                         64000000.0]
        [ 1.0, 1.0,                          -480032.0]
        [ 2.0, 1.0,                           480016.0]
        [ 1.0, 2.0,                           960016.0]
        [ 2.0, 2.0,                          -960000.0]
        [ 1.0, 3.0,                          -640000.0]
        [ 2.0, 3.0,                           640000.0]];

degree = [   0,     0;
             2,     0;
             0,     2;
             1,     1;
             1,     0;
             0,     1];
         
constraints = cell(num,1);

% This is the global optimum of El-Attar-Vidyasagar-Dutta
% f(x) = 1.712780354
feasible_point = [0.51704593415;
                  0.4891428348];

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