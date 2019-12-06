function [cost,constraints] = setup_problem_matrix_P3()
cost = [0,   0,   -10;
        1,   0,    20];

g1 = [0,     0,   110;
     1,     0,  -400;
     2,     0,   400;
     0,     1,   -20];
  
g2 = [[   0,   0,                              119]
        [ 1.0,   0,                           -680]
        [ 2.0,   0,                           1280]
        [ 3.0,   0,                           -800]
        [   0, 1.0,                              2]];
  
constraints = {g1;g2};

end

