#include "mex.h"
#include "poly2BC.h"
#include "BCDATA.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double* opt_data = (double*)mxGetData(prhs[0]);
	uint32_t opt_m = (uint32_t)mxGetM(prhs[0]);
	uint32_t opt_n = (uint32_t)mxGetN(prhs[0]);

	poly opt(opt_m - 1, opt_n, opt_data);

	double* con_data = (double*)mxGetData(prhs[1]);
	uint32_t con_m = (uint32_t)mxGetM(prhs[1]);
	uint32_t con_n = (uint32_t)mxGetN(prhs[1]);

	double* con_len = (double*)mxGetData(prhs[2]);
	uint32_t conNum = (uint32_t)mxGetM(prhs[2]);
	double* con_data_shift = con_data;

	if (conNum == 1 && con_len[0] == 0) {
		conNum = 0;
	}

	double* equ_data = (double*)mxGetData(prhs[3]);
	uint32_t equ_m = (uint32_t)mxGetM(prhs[3]);
	uint32_t equ_n = (uint32_t)mxGetN(prhs[3]);

	double* equ_len = (double*)mxGetData(prhs[4]);
	uint32_t equNum = (uint32_t)mxGetM(prhs[4]);
	double* equ_data_shift = equ_data;

	uint32_t options = (uint32_t)(*(double*)mxGetData(prhs[5]));
	bool verboseMode = options & 1;
	bool memoryRecordMode = (options >> 1) & 1;

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
	int ba_exitflag = b.solve(false, verboseMode, memoryRecordMode);

	if (memoryRecordMode) {
		nlhs = 3;
	}
	else {
		nlhs = 2;
	}

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

	plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	double* output_2 = (double*)mxGetData(plhs[1]);
	output_2[0] = (double)b.target_accuracy;

	if (memoryRecordMode) {
		plhs[2] = mxCreateNumericMatrix(2 * MAX_ITER_NUM * b.numDimension, 1, mxUINT32_CLASS, mxREAL);
		uint32_t* output_3 = (uint32_t*)mxGetData(plhs[2]);
		for (uint32_t i = 0; i < 2 * MAX_ITER_NUM * b.numDimension; i++) {
			output_3[i] = (double)b.numUnit_array[i];
		}
	}
}