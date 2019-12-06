function [cost,constraints] = setup_problem_matrix_P1(num)
cost = [1,     0,    -3;
        0,     1,    -4];

g1 = [0,     0,    -2;
     2,     0,   -72;
     3,     0,   216;
     4,     0,  -162;
     0,     1,     4];
  
g2 = [0,     0,   -36;
     1,     0,   288;
     2,     0,  -792;
     3,     0,   864;
     4,     0,  -324;
     0,     1,     4];
  
constraints = {g1;g2};

end

