function [cost,constraints] = setup_problem_matrix_P4(num)
cost = [1     0     0    -4;
        0     1     0    10;
        0     0     1    -3];

g1 = [0     0     0    -4;
      1     0     0     2;
      0     1     0    10;
      0     0     1     3];
  
g2 = [0     0     0    -6;
      0     1     0    30;
      0     0     1     3];
  
g3 = [0     0     0   -24;
     1     0     0    40;
     2     0     0   -16;
     0     1     0   -90;
     0     2     0  -200;
     1     1     0    80;
     0     0     1    39;
     0     0     2   -18;
     1     0     1   -24;
     0     1     1    60];
  
constraints = {g1;g2;g3};

end

