% setup your problem and estimate the memory usage here
% make sure that it doesn't exceed your gpu memory

clear; clc;

numDimension = 3;

% number of pre-allocated patches
twod_MAX_UNIT_NUM = 5000;
threed_MAX_UNIT_NUM = 20000;
fourd_MAX_UNIT_NUM = 100000;

if (numDimension == 2) 
    MAX_UNIT_NUM = twod_MAX_UNIT_NUM;
else
    if (numDimension == 3)
        MAX_UNIT_NUM = threed_MAX_UNIT_NUM;
    else 
        MAX_UNIT_NUM = fourd_MAX_UNIT_NUM;
    end
end

% information of polynomials
cost_degree = [2,2,2];

numCons = 128;
constraint_degree = [8,8,8];

numEqus = 0;
equality_constraint_degree = [];

opt_unitLength = prod(cost_degree + 1);
con_unitLength = prod(constraint_degree + 1);
equ_unitLength = prod(equality_constraint_degree + 1);

% compute memory usage
float_size = 4;
uint32_size = 4;
bool_size = 1/8;
char_size = 1;

mem_opt_BC = MAX_UNIT_NUM * opt_unitLength * float_size;
mem_pd_BC = numDimension * MAX_UNIT_NUM * opt_unitLength * float_size;
mem_pdValue = numDimension * MAX_UNIT_NUM * float_size;
mem_con_BC = numCons * MAX_UNIT_NUM * con_unitLength * float_size;
mem_equ_BC = numEqus * MAX_UNIT_NUM * equ_unitLength * float_size;
mem_interval = MAX_UNIT_NUM * numDimension * uint32_size;
mem_bdMin = MAX_UNIT_NUM * float_size;
mem_bdMax = MAX_UNIT_NUM * float_size;
mem_pdFlag = numDimension * MAX_UNIT_NUM * bool_size;
mem_dFlag = MAX_UNIT_NUM * bool_size;
mem_consFlag = numCons * MAX_UNIT_NUM * char_size;
mem_intFlag = MAX_UNIT_NUM * char_size;
mem_equsFlag = numEqus * MAX_UNIT_NUM * bool_size;
mem_eFlag = MAX_UNIT_NUM * bool_size;
mem_elimPos = MAX_UNIT_NUM * uint32_size;
mem_savePos = MAX_UNIT_NUM * uint32_size;

total_mem_usage = mem_opt_BC + mem_pd_BC + mem_pdValue + mem_con_BC + mem_equ_BC + mem_interval + mem_bdMin + mem_bdMax + mem_pdFlag + mem_dFlag + mem_consFlag + mem_intFlag + mem_equsFlag + mem_eFlag + mem_elimPos + mem_savePos;

fprintf("Total memory usage in GB: %f\n", total_mem_usage / 1024^3);


