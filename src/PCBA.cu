#include "mex.h"
#include "poly2BC.h"
#include "BCDATA.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double* opt_data = mxGetPr(prhs[0]);
	uint16_t opt_m = (uint16_t)mxGetM(prhs[0]);
	uint8_t opt_n = (uint8_t)mxGetN(prhs[0]);

	poly opt(opt_m - 1, opt_n, opt_data);

	double* con_data = mxGetPr(prhs[1]);
	uint16_t con_m = (uint16_t)mxGetM(prhs[1]);
	uint8_t con_n = (uint8_t)mxGetN(prhs[1]);

	double* con_len = mxGetPr(prhs[2]);
	uint16_t conNum = (uint16_t)mxGetM(prhs[2]);
	double* con_data_shift = con_data;

	if (conNum == 1 && con_len[0] == 0) {
		conNum = 0;
	}

	double* equ_data = mxGetPr(prhs[3]);
	uint16_t equ_m = (uint16_t)mxGetM(prhs[3]);
	uint8_t equ_n = (uint8_t)mxGetN(prhs[3]);

	double* equ_len = mxGetPr(prhs[4]);
	uint16_t equNum = (uint16_t)mxGetM(prhs[4]);
	double* equ_data_shift = equ_data;

	if (equNum == 1 && equ_len[0] == 0) {
		equNum = 0;
	}

	std::vector<poly> cons;
	cons.reserve(conNum);

	for (int i = 0; i < conNum; i++) {
		cons.emplace_back(con_m - 1, con_len[i], con_data_shift);
		con_data_shift += (int)con_len[i] * con_m;
	}

	std::vector<poly> equs;
	equs.reserve(equNum);

	for (int i = 0; i < equNum; i++) {
		equs.emplace_back(equ_m - 1, equ_len[i], equ_data_shift);
		equ_data_shift += (int)equ_len[i] * equ_m;
	}

	BC b(&opt, conNum, cons.data(), equNum, equs.data());
	int ba_exitflag = b.solve(false, false);

	nlhs = 3;
	if (ba_exitflag != -12345)
	{
		plhs[0] = mxCreateNumericMatrix((mwSize)b.numDimension, 1, mxDOUBLE_CLASS, mxREAL);
		double* output = (double*)mxGetData(plhs[0]);
		for (int dim = 0; dim < b.numDimension; dim++) {
			output[dim] = (double)(b.final_result[dim]);
		}
	}
	else
	{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double* output = (double*)mxGetData(plhs[0]);
		output[0] = (double)ba_exitflag;
	}

	plhs[1] = mxCreateNumericMatrix(2 * MAX_ITER_NUM * b.numDimension, 1, mxUINT32_CLASS, mxREAL);
	uint32_t* output_2 = (uint32_t*)mxGetData(plhs[1]);
	for(uint32_t i = 0; i < 2 * MAX_ITER_NUM * b.numDimension; i++){
		output_2[i] = (double)b.numUnit_array[i];
	}

	plhs[2] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	double* output_3 = (double*)mxGetData(plhs[2]);
	output_3[0] = (double)b.target_accuracy;
}