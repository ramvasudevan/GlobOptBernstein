/*
Constrained Global Optimization Using Bernstein Algorithm

Copyright: ROAHMLAB, University of Michigan

Arthor: Bohao Zhang (jimzhang@umich.edu)

Insturction: This file defines a Bernstein Algorithm solver class and all the cuda kernel functions that it needs.

Basic Introduction to the data structure and its motivation:
To make the data as 'parallelized' as possible, The Bernstein coefficients are not stored in any order.
As a result, additional to the data sturture storing the Bernstein coefficients, we need another array to store the 'position' of this patch
inside the variable space, so that we are able to compute the solution corresponding with a specific patch.
For the Bernstein algorithm, the 'position' is represented by a box \{(x_1,...,x_m) \in R^m | 0<=a_i<=x_i<=b_i<=1, i=1,...,m\}
This 'position' array will, hence, be performed with exactly the same operation together with the data sturture storing the Bernstein coefficients.
In general, we define a unit to be a union of the following elements:
	 1. a set of integers indicating the position (box) in variable space
	 2. Bernstein coefficients of the problem polynomial over that box
	 3. Bernstein coefficients of the inequality constraint polynomials over that box
	 4. Bernstein coefficients of the equality constraint polynomials over that box
	 5. lower and upper bound of the problem polynomial over that box (computed from 2.)
	 6. Relevant flags labelling the feasibility of that box (computed from 3. and 4.)
All the information is promised to be stored in the same index in no matter which data structure.
*/

#ifndef _BCDATA_H_
#define _BCDATA_H_

#include "poly2BC.h"

// maximum of iteration number, should NOT be larger than 32 since position info is stored in uint32_t
#define MAX_ITER_NUM 30

// used for shared memory pre-allocation, the maximum size of a Bernstein patch
#define MAX_UNIT_LENGTH 1000

// used for global memory pre-allocation, the maximum number of Bernstein patches for different dimension problems
#define twod_MAX_UNIT_NUM  5000
#define threed_MAX_UNIT_NUM 40000
#define fourd_MAX_UNIT_NUM 100000

#define STOPPING_CRITERIA 1e-7

class BC {
public:
	/*
	REQUIRES:
		opt_in     -> the polynomial needed to be optimized
		numCons_in -> number of inequality constraints. has to be smaller than 256 here. change data type if needed
		cons_in    -> an array of inequality constraint polynomials
		numEqus_in -> number of equality constraints. has to be smaller than 256 here. change data type if needed
		equs_in    -> an array of equality constraint polynomials
	EFFECT: defined constructor to initialize the optimizer
	*/
	BC(poly* opt_in, uint8_t numCons_in, poly* cons_in, uint8_t numEqus_in, poly* equs_in);

	/*
	EFFECT: default destructor
	*/
	~BC();

	/*
	EFFECT: print data related on the screen. used for debug
	*/
	void debug_print();

	/*
	REQUIRES:
		dim -> the index of the dimension where subdivision is performed. dim < numDimension
	EFFECT: divide the patches (problem polynomial and constraint polynomials) into two equally. 
			
			[info for 1 st unit,         info for 2 nd unit,         ..., info for n th unit]
													|
													|
											 after dilation
													|
													V
		    [info for left of 1 st unit, info for left of 2 nd unit, ..., info for n th left of unit, info for right of 1 st unit, info for right of 2 nd unit, ..., info for n th right of unit]
	*/
	void dilation(uint8_t dim);

	/*
	EFFECT: find the lower and upper bound of the patches (problem polynomial and constraint polynomials).
			Thus we are able to judge the feasibility of the patch
	*/
	void findFlag();

	/*
	EFFECT: locate the lowest feasible upper bound as estimated minimum.
		    eliminate the infeasible patches and those whose lower bound is larger than the estimated minimum
			To be more exact, the patches needed to be eliminated will be replaced by the patches needed to be saved

			Example: 
			[info for 1 st unit, info for 2 nd unit, info for 3 rd unit, info for 4 th unit, info for 5 th unit, info for 6 th unit, info for 7 th unit]
			 eliminate 1st       save 1st            eliminate 2nd       save 2nd            save 3rd            eliminate 3rd       save 4th
																					|
																					|
																			 after elimination
																					|
																					V
			[info for 7 th unit, info for 2 nd unit, info for 5 th unit, info for 4 th unit]	
			 save 4th            save 1st            save 3rd            save 2nd

			The patches needed to be saved remain at their original position. And The patches need to be eliminated will be replaced by saved patches.
			In that way, elimination can be performed parallizedly. And that's why we need another array to store the interval where the patches lie in.

	*/
	void eliminate();

	/*
	EFFECT: All the iterations have been completed. Convert that unit into the result we need (the optimum and the box relevant)
	*/
	void finalResult();

	/*
	REQUIRES:
		debugMode   -> turn on to print the data on the screen at each iteration. 
					   Note that you have to modify the constructor otherwise memory leak would take place!!!!!!!
					   default: turn off
		verboseMode -> print number of patches left at each iteration on the screen
					   default: turn on
	EFFECT: the optimizer to call. put everything together as described in the paper. store the result into final_result
	RETURN:
		1      -> optimization successfully finished
		-12345 -> the problem is infeasible
		-54321 -> too much memory is needed. the optimizer didn't compplete all the iterations required in accuracy. A less precise result is given
	*/
	int solve(bool debugMode = false, bool verboseMode = true);

	/*
	The maximum number of Benrstein units possible in the program. used for pre-allocating memory.
	This number shall vary according to the dimension of the problem since the number of existed patches may increase exponentially with the dimension.
	In this project, we set
		2d: 5000
		3d: 40000
		4d: 100000
	This limit is proved to work for most of the optimization problems.
	*/
	uint32_t MAX_UNIT_NUM;

	// a pointer to the problem polynomial needed to be optimized
	poly* opt;

	// an array of pointers to the inequality constraint polynomials
	poly* cons;
	uint8_t numCons;

	// an array of pointers to the equality constraint polynomials
	poly* equs;
	uint8_t numEqus;

	// number of dimension of the problem
	uint8_t numDimension;

	// an array of (maximum) degrees of x_i for (problem, ineq con, eq con) polynomials
	uint8_t* opt_degree;
	uint8_t* con_degree;
	uint8_t* equ_degree;

	// the length of a Bernstein unit. In other word, the number of Bernstein coefficients in a unit.
	uint16_t opt_unitLength;
	uint16_t con_unitLength;
	uint16_t equ_unitLength;

	/* 
	the array of Bernstein coefficients of the problem polynomial. 
	Recall that Bernstein coefficients actually are represented by the polynomial estimation over a grid of a m-dimensional cube, where m is the dimension of the problem.
	As a result, we will need to store a m-dimensional cube into a 1-dimensional array. The structure of the array is like this
	[                                                                                                                                                                                                                                                                          ]
	 <---- 1st dim ----> ... <---- 1st dim ----> ... <---- 1st dim ----> ... <---- 1st dim ----> ... <---- 1st dim ----> ... <---- 1st dim ----> ... <---- 1st dim ----> ... <---- 1st dim ---->
	 <----          1st * 2nd dim          ----> ... <----          1st * 2nd dim          ----> ... <----          1st * 2nd dim          ----> ... <----          1st * 2nd dim          ---->
	 <----                              1st * 2nd * 3rd dim                                ----> ... <----                              1st * 2nd * 3rd dim                                ---->
	 .
	 .
	 .                                                                                                                                                                                   ......
	 <----                                                                                      1st unit                                                                                   ---->, <---- 2nd unit ---->, ..., <---- nth unit ---->                                  
	 <----                                                                                   opt_unitLength                                                                                ---->
	*/
	float* opt_BC;
	float* dev_opt_BC;

	/*
	the array of Bernstein coefficients of the partial derivatives of the problem polynomial. The structure of the array is like this
	[                                                                                                                                                                               ]
	 <---- Bernstein coefs of df / dx_1 ----> ... <---- Bernstein coefs of df / dx_m ----> ... <---- Bernstein coefs of df / dx_1 ----> ... <---- Bernstein coefs of df / dx_m ---->
	 <----        opt_unitLength        ---->
	 <----                                  1st unit                                 ----> ... <----                                  nth unit                                 ---->    
	*/
	float* pd_BC;
	float* dev_pd_BC;

	/*
	the array of Bernstein coefficients of the inequality constraint polynomial. The structure of the array is like this
	[                                                                                                                                                                           ]
	 <---- Bernstein coefs of 1st cons ----> ... <---- Bernstein coefs of nth cons ----> ... <---- Bernstein coefs of 1st cons ----> ... <---- Bernstein coefs of nth cons ---->
	 <----       con_unitLength        ---->
	 <----                                1st unit                                 ----> ... <----                                 nth unit                                ---->
	*/
	float* con_BC;
	float* dev_con_BC;

	/*
	the array of Bernstein coefficients of the equality constraint polynomial. The structure of the array is like this
	[                                                                                                                                                                           ]
	 <---- Bernstein coefs of 1st cons ----> ... <---- Bernstein coefs of nth cons ----> ... <---- Bernstein coefs of 1st cons ----> ... <---- Bernstein coefs of nth cons ---->
	 <----       equ_unitLength        ---->
	 <----                                1st unit                                 ----> ... <----                                 nth unit                                ---->
	*/
	float* equ_BC;
	float* dev_equ_BC;

	/*
	the array indicating the corresponding box of that unit. 
	It can be represented by m integers where m is the number of dimension, one number for one dimension.
	The exact numerical value of the box can be computed as [x / (2 ^ (number of iteration on this dim)), (x + 1) / (2 ^ (number of iteration on this dim))]
	All the integers are initialized to 0 before the solver starts.

	Example:
	When two iterations have been completed and the current dimension is x_2 , one box becomes [3, 0,2].
	It means the box [0.375, 0.5] * [0, 0.125] * [0.5, 0.75] in a 3-dimension space.
	This is because x_1 and x_2 has been divided for three times while x_3 has only been divided for two times.

	Such representation allows us to compute subdivision quickly.
	x -> original integer on this dimension
	(x << 1) -> left sub-box integer on this dimension
	(x << 1) + 1 -> right sub-box integer on this dimension
	*/
	uint32_t* interval;
	uint32_t* dev_interval;

	/*
	As the iteration performed on different dimension could be different if the program breaks too early,
	we will also need one array recording the number of iterations having been performed on each dimension,
	so that we could compute the solutions at the end of the program.
	*/
	uint8_t* int_iter;

	/*
	The constant coefficient of each partial derivatives. 
	will be needed to judge whether the bound of the partial derivative includes 0.
	generated from poly::partialDerivative
	*/
	float* pdValue;
	float* dev_pdValue;

	/*
	A floating number array storing the exact results of boxes represented by integer numbers.
	This will be needed for the last step of the program, after some stopping criteria has been satisfied.
	*/
	float* intervalRes;
	float* dev_intervalRes;

	/*
	number of patches remained
	We record the apex of the number of patches during the optimization to observe the performance of the algorithm in term of memory usage
	*/
	uint32_t numUnit;
	uint32_t apex_numUnit;

	/*
	the lower bound of a Bernstein patch
	in other words, the minimum of Bernstein coefficients over a certain box
	used for cut-off test
	*/
	float* bdMin;
	float* dev_bdMin;

	/*
	the upper bound of a Bernstein patch
	in other words, the maximum of Bernstein coefficients over a certain box
	used for cut-off test
	*/
	float* bdMax;
	float* dev_bdMax;

	/*
	a flag labelling whether a Bernstein patch for partial derivative includes 0 along a certain dimension
	in other words, lower bound <= 0 <= upper bound
	Using this flag, we can tell whether there's a local extrema along that dimension. 
	*/
	bool* pdFlag;
	bool* dev_pdFlag;

	/*
	a flag labelling whether there's local extremas along all the dimension. 
	Using this flag, we can tell whether the function is strictly increasing or decreasing along all dimension.
	If not, and the box is not on the boundary of any constriants, we will eliminate this box since it will definitely not include a global optimum.
	*/
	bool* dFlag;
	bool* dev_dFlag;

	/*
	a flag labelling whether the box satisfies a certain inequality constraint
	lower bound >= 0               --> consFlag = 0 --> not satisfy
	lower bound < 0 < upper bound  --> consFlag = 1 --> not strictly satisfy
	upper bound <= 0               --> consFlag = 2 --> strictly satisfy
	*/
	char* consFlag;
	char* dev_consFlag;

	/*
	a flag labelling whether the box satisfies all the inequality constraints
	in other words, and (&) operation over all consFlags on that box
	*/
	char* intFlag;
	char* dev_intFlag;

	/*
	a flag labelling whether the box satisfies a certain equality constraint
	lower bound > 0 || upper bound < 0 --> equsFlag = 0 --> not satisfy
	lower bound <= 0 <= upper bound    --> equsFlag = 1 --> satisfy
	*/
	bool* equsFlag;
	bool* dev_equsFlag;

	/*
	a flag labelling whether the box satisfies all the equality constraints
	in other words, and (&) operation over all equsFlag on that box
	*/
	bool* eFlag;
	bool* dev_eFlag;

	/*
	the index of boxes needed to be eliminated
	*/
	uint32_t* elimPos;
	uint32_t* dev_elimPos;

	/*
	the index of boxes needed to be saved
	*/
	uint32_t* savePos;
	uint32_t* dev_savePos;

	/*
	current estimated optimum
	*/
	float estiMin;

	/*
	current iteration and dimension
	*/
	uint8_t iter;
	uint8_t dim;

	/*
	the optimal polynomial values over the center certain boxes
	After all the iterations finished, there will be still some boxes remained, satisfying the stopping criteria.
	We compute the exact value of the optimal polynomial over the center of them and then find the minimum among them,
	so that the results returned will be more accurate.
	*/
	float* candidates;

	/*
	the array storing the final results, i.e. the minimum from candidates
	and its corresponding index in the array
	*/
	float* final_result;
	uint32_t final_index;

	/*
	the target accuracy = STOPPING_CRITERIA * (upper bound of the initial Bernstein patch - lower bound of the initial Bernstein patch)
	Our paper includes more details about this.
	*/
	float target_accuracy;

	/*
	the accuracy of the current iteration = minimum of all the upper bounds - minimum of all the lower bounds
	Our paper includes more details about this.
	*/
	float estimated_accuracy;

	// a flag for whether the stopping criteria has been satisfied, so that the loop can end
	bool last;

	float* debug;
	float* dev_debug;
};

__global__ void biBCdilationKernelPart1Forx1(float* target_BC, uint8_t x1degree, uint16_t unitLength);

__global__ void biBCdilationKernelPart1Forx2(float* target_BC, uint8_t x2degree, uint16_t unitLength);

__global__ void triBCdilationKernelPart1Forx1(float* target_BC, uint8_t x1degree, uint16_t unitLength);

__global__ void triBCdilationKernelPart1Forx2(float* target_BC, uint8_t x2degree, uint16_t unitLength);

__global__ void triBCdilationKernelPart1Forx3(float* target_BC, uint8_t x2degree, uint16_t unitLength);

__global__ void quadBCdilationKernelPart1Forx1(float* target_BC, uint8_t x1degree, uint16_t unitLength);

__global__ void quadBCdilationKernelPart1Forx2(float* target_BC, uint8_t x2degree, uint16_t unitLength);

__global__ void quadBCdilationKernelPart1Forx3(float* target_BC, uint8_t x3degree, uint16_t unitLength);

__global__ void quadBCdilationKernelPart1Forx4(float* target_BC, uint8_t x4degree, uint16_t unitLength);

__global__ void BCdilationKernelPart2(uint32_t* target_interval, bool* target_pdFlag, uint8_t dim, float* target_pdValue);

__global__ void BCdilationKernelPart3(char* target_consFlag);

__global__ void BCfindFlagKernel(char* conFlag, float* BC);

__global__ void BCfindIntFlagKernel(char* intFlag, char* conFlag, uint8_t numCons);

__global__ void BCfindEquFlagKernel(bool* equFlag, float* BC);

__global__ void BCfindEFlagKernel(bool* eFlag, bool* equsFlag, uint8_t numEqus);

__global__ void BCfindBoundKernel(float* bdMin, float* bdMax, float* BC);

__global__ void BCfindDerivativeKernel(bool* pdFlag, float* BC, float* pdValue);

__global__ void biBCfinddFlagKernel(bool* dFlag, bool* pdFlag, uint8_t numDimension, uint32_t* interval, uint8_t iter_0, uint8_t iter_1);

__global__ void triBCfinddFlagKernel(bool* dFlag, bool* pdFlag, uint8_t numDimension, uint32_t* interval, uint8_t iter_0, uint8_t iter_1, uint8_t iter_2);

__global__ void quadBCfinddFlagKernel(bool* dFlag, bool* pdFlag, uint8_t numDimension, uint32_t* interval, uint8_t iter_0, uint8_t iter_1, uint8_t iter_2, uint8_t iter_3);

__global__ void BCeliminateKernelPart1(float* target_BC, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum);

__global__ void BCeliminateKernelPart2(uint32_t* target_interval, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum);

__global__ void BCeliminateKernelPart3(float* target_BC, char* target_consFlag, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum);

__global__ void BCeliminateKernelPart4(float* target_BC, bool* target_pdFlag, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum);

__global__ void biBCfinalResultKernel(float* target_intervalRes, uint32_t* interval, uint8_t iter_1, uint8_t iter_2);

__global__ void triBCfinalResultKernel(float* target_intervalRes, uint32_t* interval, uint8_t iter_1, uint8_t iter_2, uint8_t iter_3);

__global__ void quadBCfinalResultKernel(float* target_intervalRes, uint32_t* interval, uint8_t iter_1, uint8_t iter_2, uint8_t iter_3, uint8_t iter_4);

#endif
