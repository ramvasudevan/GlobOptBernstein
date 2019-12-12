#ifndef _BCDATA_CPP_
#define _BCDATA_CPP_

#include "BCDATA.h"

BC::BC(poly* opt_in, uint8_t numCons_in, poly* cons_in, uint8_t numEqus_in, poly* equs_in) {
	opt = opt_in;
	cons = cons_in;
	numCons = numCons_in;
	equs = equs_in;
	numEqus = numEqus_in;
	numUnit = 1;
	apex_numUnit = 1;
	numDimension = opt_in->numDimension;

	if (numDimension == 2) {
		MAX_UNIT_NUM = twod_MAX_UNIT_NUM;
	}
	else if (numDimension == 3) {
		MAX_UNIT_NUM = threed_MAX_UNIT_NUM;
	}
	else {
		MAX_UNIT_NUM = fourd_MAX_UNIT_NUM;
	}

	opt_degree = new uint8_t[numDimension];
	if (numCons > 0) con_degree = new uint8_t[numDimension];
	else con_degree = nullptr;
	if (numEqus > 0) equ_degree = new uint8_t[numDimension];
	else equ_degree = nullptr;

	opt_unitLength = con_unitLength = equ_unitLength = 1;
	for (uint8_t i = 0; i < numDimension; i++)
	{
		opt_degree[i] = opt_in->maxDegree[i];
		opt_unitLength *= opt_degree[i];

		if (numCons > 0) {
			con_degree[i] = 0;
			for (uint8_t j = 0; j < numCons; j++) {
				if (cons_in[j].maxDegree[i] > con_degree[i]) {
					con_degree[i] = cons_in[j].maxDegree[i];
				}
			}
			con_unitLength *= con_degree[i];
		}

		if (numEqus > 0) {
			equ_degree[i] = 0;
			for (uint8_t j = 0; j < numEqus; j++) {
				if (equs_in[j].maxDegree[i] > equ_degree[i]) {
					equ_degree[i] = equs_in[j].maxDegree[i];
				}
			}
			equ_unitLength *= equ_degree[i];
		}
	}

	//opt_BC = new float[MAX_UNIT_NUM * opt_unitLength];
	opt_BC = new float[opt_unitLength];
	memset(opt_BC, 0, opt_unitLength * sizeof(float));

	uint16_t index, currentPos, currentDegree;
	float BCterm;

	for (uint16_t i = 0; i < opt_in->numTerms; i++) {
		for (uint16_t j = 0; j < opt_unitLength; j++) {
			index = j;
			BCterm = 1;
			for (int k = 0; k < numDimension; k++) {
				currentDegree = index % opt_degree[k];
				currentPos = i * numDimension + k;
				BCterm *= choosenk[currentDegree][opt_in->degree[currentPos]] / choosenk[opt_degree[k] - 1][opt_in->degree[currentPos]];
				index = (index - currentDegree) / opt_degree[k];
			}
			opt_BC[j] += opt_in->coeff[i] * BCterm;
		}
	}

	cudaMalloc((void**)&dev_opt_BC, MAX_UNIT_NUM * opt_unitLength * sizeof(float));
	cudaMemcpy(dev_opt_BC, opt_BC, opt_unitLength * sizeof(float), cudaMemcpyHostToDevice);

	pdValue = new float[numDimension];
	//pd_BC = new float[MAX_UNIT_NUM * numDimension * opt_unitLength];
	pd_BC = new float[numDimension * opt_unitLength];
	memset(pd_BC, 0, numDimension * opt_unitLength * sizeof(float));

	uint16_t pdOffset = 0;
	poly* pd = nullptr;
	float pdValueBuf = 0;
	for (uint8_t pdID = 0; pdID < numDimension; pdID++) {
		opt->partialDerivative(pd, pdValueBuf, pdID);
		for (uint16_t i = 0; i < pd->numTerms; i++) {
			for (uint16_t j = 0; j < opt_unitLength; j++) {
				index = j;
				BCterm = 1;
				for (int k = 0; k < numDimension; k++) {
					currentDegree = index % opt_degree[k];
					currentPos = i * numDimension + k;
					BCterm *= choosenk[currentDegree][pd->degree[currentPos]] / choosenk[opt_degree[k] - 1][pd->degree[currentPos]];
					index = (index - currentDegree) / opt_degree[k];
				}
				pd_BC[pdOffset + j] += pd->coeff[i] * BCterm;
			}
		}
		pdOffset += opt_unitLength;
		delete pd;

		pdValue[pdID] = pdValueBuf;
	}

	cudaMalloc((void**)&dev_pd_BC, numDimension * MAX_UNIT_NUM * opt_unitLength * sizeof(float));
	cudaMemcpy(dev_pd_BC, pd_BC, numDimension * opt_unitLength * sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&dev_pdValue, numDimension * MAX_UNIT_NUM * sizeof(float));
	cudaMemcpy(dev_pdValue, pdValue, numDimension * sizeof(float), cudaMemcpyHostToDevice);

	if (numCons > 0) {
		//con_BC = new float[numCons * MAX_UNIT_NUM * con_unitLength];
		con_BC = new float[numCons * con_unitLength];
		memset(con_BC, 0, numCons * con_unitLength * sizeof(float));

		uint16_t conOffset = 0;

		for (uint8_t conID = 0; conID < numCons; conID++) {
			for (uint16_t i = 0; i < cons[conID].numTerms; i++) {
				for (uint16_t j = 0; j < con_unitLength; j++) {
					index = j;
					BCterm = 1;
					for (int k = 0; k < numDimension; k++) {
						currentDegree = index % con_degree[k];
						currentPos = i * numDimension + k;
						BCterm *= choosenk[currentDegree][cons[conID].degree[currentPos]] / choosenk[con_degree[k] - 1][cons[conID].degree[currentPos]];
						index = (index - currentDegree) / con_degree[k];
					}
					con_BC[conOffset + j] += cons[conID].coeff[i] * BCterm;
				}
			}
			conOffset += con_unitLength;
		}

		cudaMalloc((void**)&dev_con_BC, numCons * MAX_UNIT_NUM * con_unitLength * sizeof(float));
		cudaMemcpy(dev_con_BC, con_BC, numCons * con_unitLength * sizeof(float), cudaMemcpyHostToDevice);
	}
	else {
		con_BC = nullptr;
		dev_con_BC = nullptr;
	}

	if (numEqus > 0) {
		//equ_BC = new float[numEqus * MAX_UNIT_NUM * equ_unitLength];
		equ_BC = new float[numEqus * equ_unitLength];
		memset(equ_BC, 0, numEqus * equ_unitLength * sizeof(float));

		uint16_t equOffset = 0;

		for (uint8_t equID = 0; equID < numEqus; equID++) {
			for (uint16_t i = 0; i < equs[equID].numTerms; i++) {
				for (uint16_t j = 0; j < equ_unitLength; j++) {
					index = j;
					BCterm = 1;
					for (int k = 0; k < numDimension; k++) {
						currentDegree = index % equ_degree[k];
						currentPos = i * numDimension + k;
						BCterm *= choosenk[currentDegree][equs[equID].degree[currentPos]] / choosenk[equ_degree[k] - 1][equs[equID].degree[currentPos]];
						index = (index - currentDegree) / equ_degree[k];
					}
					equ_BC[equOffset + j] += equs[equID].coeff[i] * BCterm;
				}
			}
			equOffset += equ_unitLength;
		}

		cudaMalloc((void**)&dev_equ_BC, numEqus * MAX_UNIT_NUM * equ_unitLength * sizeof(float));
		cudaMemcpy(dev_equ_BC, equ_BC, numEqus * equ_unitLength * sizeof(float), cudaMemcpyHostToDevice);
	}
	else {
		equ_BC = nullptr;
		dev_equ_BC = nullptr;
	}

	//interval = new uint32_t[MAX_UNIT_NUM * numDimension];
	interval = new uint32_t[numDimension];
	memset(interval, 0, numDimension * sizeof(uint32_t));
	cudaMalloc((void**)&dev_interval, MAX_UNIT_NUM * numDimension * sizeof(uint32_t));
	cudaMemcpy(dev_interval, interval, numDimension * sizeof(uint32_t), cudaMemcpyHostToDevice);

	int_iter = new uint8_t[numDimension];
	memset(int_iter, 0, numDimension * sizeof(uint8_t));

	bdMin = new float[MAX_UNIT_NUM];
	cudaMalloc((void**)&dev_bdMin, MAX_UNIT_NUM * sizeof(float));

	bdMax = new float[MAX_UNIT_NUM];
	cudaMalloc((void**)&dev_bdMax, MAX_UNIT_NUM * sizeof(float));

	//pdFlag = new bool[MAX_UNIT_NUM * numDimension];
	pdFlag = new bool[numDimension];
	for (uint8_t i = 0; i < numDimension; i++) {
		pdFlag[i] = true;
	}
	cudaMalloc((void**)&dev_pdFlag, numDimension * MAX_UNIT_NUM * sizeof(bool));
	cudaMemcpy(dev_pdFlag, pdFlag, numDimension * sizeof(bool), cudaMemcpyHostToDevice);

	dFlag = new bool[MAX_UNIT_NUM];
	cudaMalloc((void**)&dev_dFlag, MAX_UNIT_NUM * sizeof(bool));

	if (numCons > 0) {
		//consFlag = new char[numCons * MAX_UNIT_NUM];
		consFlag = new char[numCons];
		for (uint8_t i = 0; i < numCons; i++) {
			consFlag[i] = 1;
		}
		cudaMalloc((void**)&dev_consFlag, numCons * MAX_UNIT_NUM * sizeof(char));
		cudaMemcpy(dev_consFlag, consFlag, numCons * sizeof(char), cudaMemcpyHostToDevice);

		intFlag = new char[MAX_UNIT_NUM];
		cudaMalloc((void**)&dev_intFlag, MAX_UNIT_NUM * sizeof(char));
	}
	else {
		consFlag = nullptr;
		dev_consFlag = nullptr;
		intFlag = nullptr;
		dev_intFlag = nullptr;
	}

	if (numEqus > 0) {
		//equsFlag = new bool[numEqus * MAX_UNIT_NUM];
		equsFlag = new bool[numEqus];
		cudaMalloc((void**)&dev_equsFlag, numEqus * MAX_UNIT_NUM * sizeof(bool));

		eFlag = new bool[MAX_UNIT_NUM];
		cudaMalloc((void**)&dev_eFlag, MAX_UNIT_NUM * sizeof(bool));
	}
	else {
		equsFlag = nullptr;
		dev_equsFlag = nullptr;
		eFlag = nullptr;
		dev_eFlag = nullptr;
	}


	elimPos = new uint32_t[MAX_UNIT_NUM];
	cudaMalloc((void**)&dev_elimPos, MAX_UNIT_NUM * sizeof(uint32_t));

	savePos = new uint32_t[MAX_UNIT_NUM];
	cudaMalloc((void**)&dev_savePos, MAX_UNIT_NUM * sizeof(uint32_t));

	intervalRes = nullptr;
	dev_intervalRes = nullptr;

	candidates = nullptr;

	//debug = new float[MAX_UNIT_NUM * numCons * 2];
	//cudaMalloc((void**)&dev_debug, MAX_UNIT_NUM * numCons * 2 * sizeof(float));

	final_index = 0;

	final_result = nullptr;
}

BC::~BC() {
	delete[] opt_degree;

	if (numCons > 0) delete[] con_degree;

	if (numEqus > 0) delete[] equ_degree;

	delete[] opt_BC;
	cudaFree(dev_opt_BC);

	delete[] pd_BC;
	cudaFree(dev_pd_BC);

	delete[] pdValue;
	cudaFree(dev_pdValue);

	if (numCons > 0) {
		delete[] con_BC;
		cudaFree(dev_con_BC);
	}

	if (numEqus > 0) {
		delete[] equ_BC;
		cudaFree(dev_equ_BC);
	}

	delete[] interval;
	cudaFree(dev_interval);

	delete[] int_iter;

	delete[] bdMin;
	cudaFree(dev_bdMin);

	delete[] bdMax;
	cudaFree(dev_bdMax);

	delete[] pdFlag;
	cudaFree(dev_pdFlag);

	delete[] dFlag;
	cudaFree(dev_dFlag);

	if (numCons > 0) {
		delete[] consFlag;
		cudaFree(dev_consFlag);

		delete[] intFlag;
		cudaFree(dev_intFlag);
	}

	if (numEqus > 0) {
		delete[] equsFlag;
		cudaFree(dev_equsFlag);

		delete[] eFlag;
		cudaFree(dev_eFlag);
	}

	delete[] elimPos;
	cudaFree(dev_elimPos);

	delete[] savePos;
	cudaFree(dev_savePos);

	if (intervalRes != nullptr) {
		delete[] intervalRes;
		cudaFree(dev_intervalRes);
	}

	if (candidates != nullptr) delete[] candidates;

	if (final_result != nullptr) delete[] final_result;

	//delete[] debug;
	//cudaFree(dev_debug);
}

void BC::debug_print()
{
	int debugNumUnit = numUnit > 10 ? 10 : numUnit;
	/*
	mexPrintf("PRINT OPT BC\n");
	cudaMemcpy(opt_BC, dev_opt_BC, debugNumUnit * opt_unitLength * sizeof(float), cudaMemcpyDeviceToHost);
	for (uint32_t i = 0; i < debugNumUnit; i++)
	{
		mexPrintf("BC %d. \n", i);
		for (int j = 0; j < opt_unitLength; j++)
		{
			mexPrintf("%.6f ", opt_BC[i * opt_unitLength + j]);
		}
		mexPrintf("\n");
	}
	mexPrintf("\n\n");
	
	mexPrintf("PRINT PD BC\n");
	cudaMemcpy(pd_BC, dev_pd_BC, debugNumUnit * numDimension * opt_unitLength * sizeof(float), cudaMemcpyDeviceToHost);
	for (uint32_t i = 0; i < debugNumUnit * numDimension; i++)
	{
		mexPrintf("BC %d. \n", i);
		for (int j = 0; j < opt_unitLength; j++)
		{
			mexPrintf("%.6f ", pd_BC[i * opt_unitLength + j]);
		}
		mexPrintf("\n");
	}
	mexPrintf("\n\n");
	
	if (numCons > 0) {
		mexPrintf("PRINT CON BC\n");
		cudaMemcpy(con_BC, dev_con_BC, debugNumUnit * numCons * con_unitLength * sizeof(float), cudaMemcpyDeviceToHost);
		for (uint32_t i = 0; i < debugNumUnit * numCons; i++)
		{
			mexPrintf("BC %d. \n", i);
			for (int j = 0; j < con_unitLength; j++)
			{
				mexPrintf("%.6f ", con_BC[i * con_unitLength + j]);
			}
			mexPrintf("\n");
		}
		mexPrintf("\n\n");
	}
	
	if (numEqus > 0) {
		mexPrintf("PRINT EQU BC\n");
		cudaMemcpy(equ_BC, dev_equ_BC, debugNumUnit * numEqus * equ_unitLength * sizeof(float), cudaMemcpyDeviceToHost);
		for (uint32_t i = 0; i < debugNumUnit * numEqus; i++)
		{
			mexPrintf("BC %d. \n", i);
			for (int j = 0; j < equ_unitLength; j++)
			{
				mexPrintf("%.6f ", equ_BC[i * equ_unitLength + j]);
			}
			mexPrintf("\n");
		}
		mexPrintf("\n\n");
	}
	*/
	if (numCons > 0) {
		cudaMemcpy(consFlag, dev_consFlag, debugNumUnit * numCons * sizeof(char), cudaMemcpyDeviceToHost);
		cudaMemcpy(intFlag, dev_intFlag, debugNumUnit * sizeof(char), cudaMemcpyDeviceToHost);
	}
	if (numEqus > 0) {
		cudaMemcpy(equsFlag, dev_equsFlag, debugNumUnit * numEqus * sizeof(bool), cudaMemcpyDeviceToHost);
		cudaMemcpy(eFlag, dev_eFlag, debugNumUnit * sizeof(bool), cudaMemcpyDeviceToHost);
	}
	cudaMemcpy(pdFlag, dev_pdFlag, debugNumUnit * numDimension * sizeof(bool), cudaMemcpyDeviceToHost);
	cudaMemcpy(dFlag, dev_dFlag, debugNumUnit * sizeof(bool), cudaMemcpyDeviceToHost);
	cudaMemcpy(interval, dev_interval, debugNumUnit * numDimension * sizeof(uint32_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(bdMin, dev_bdMin, debugNumUnit * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(bdMax, dev_bdMax, debugNumUnit * sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(debug, dev_debug, debugNumUnit * numCons * 2 * sizeof(float), cudaMemcpyDeviceToHost);

	mexPrintf("PRINT UNIT INFO\n");
	for (uint32_t i = 0; i < debugNumUnit; i++)
	{
		mexPrintf("%d. [", i);
		for (uint8_t j = 0; j < numDimension; j++) {
			mexPrintf(" %d", interval[i * numDimension + j]);
		}
		mexPrintf("]  ");

		mexPrintf("[%.6f %.6f]  | ", bdMin[i], bdMax[i]);

		if (numCons > 0) {
			for (uint8_t j = 0; j < numCons; j++) {
				mexPrintf("%d ", consFlag[i * numCons + j]);
			}
			mexPrintf(":%d | ", intFlag[i]);
		}

		if (numEqus > 0) {
			for (uint8_t j = 0; j < numEqus; j++) {
				mexPrintf("%d ", equsFlag[i * numEqus + j]);
			}
			mexPrintf(":%d | ", eFlag[i]);
		}

		for (uint8_t j = 0; j < numDimension; j++) {
			mexPrintf("%d ", pdFlag[i * numDimension + j]);
		}
		mexPrintf(":%d\n", dFlag[i]);
	}
	mexPrintf("\n");
}

void BC::dilation(uint8_t dim) {
	if (numDimension == 2) {
		if (dim == 0) {
			biBCdilationKernelPart1Forx1 << < numUnit, opt_degree[1] >> > (dev_opt_BC, opt_degree[0], opt_unitLength);
			biBCdilationKernelPart1Forx1 << < numUnit * numDimension, opt_degree[1] >> > (dev_pd_BC, opt_degree[0], opt_unitLength);
			if (numCons > 0) biBCdilationKernelPart1Forx1 << < numUnit * numCons, con_degree[1] >> > (dev_con_BC, con_degree[0], con_unitLength);
			if (numEqus > 0) biBCdilationKernelPart1Forx1 << < numUnit * numEqus, equ_degree[1] >> > (dev_equ_BC, equ_degree[0], equ_unitLength);
		}
		else {
			biBCdilationKernelPart1Forx2 << < numUnit, opt_degree[0] >> > (dev_opt_BC, opt_degree[1], opt_unitLength);
			biBCdilationKernelPart1Forx2 << < numUnit * numDimension, opt_degree[0] >> > (dev_pd_BC, opt_degree[1], opt_unitLength);
			if (numCons > 0) biBCdilationKernelPart1Forx2 << < numUnit * numCons, con_degree[0] >> > (dev_con_BC, con_degree[1], con_unitLength);
			if (numEqus > 0) biBCdilationKernelPart1Forx2 << < numUnit * numEqus, equ_degree[0] >> > (dev_equ_BC, equ_degree[1], equ_unitLength);
		}
	}
	else if (numDimension == 3) {
		if (dim == 0) {
			dim3 block1(opt_degree[1], opt_degree[2], 1);
			triBCdilationKernelPart1Forx1 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[0], opt_unitLength);
			triBCdilationKernelPart1Forx1 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[0], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[1], con_degree[2], 1);
				triBCdilationKernelPart1Forx1 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[0], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[1], equ_degree[2], 1);
				triBCdilationKernelPart1Forx1 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[0], equ_unitLength);
			}
		}
		else if (dim == 1) {
			dim3 block1(opt_degree[0], opt_degree[2], 1);
			triBCdilationKernelPart1Forx2 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[1], opt_unitLength);
			triBCdilationKernelPart1Forx2 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[1], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[0], con_degree[2], 1);
				triBCdilationKernelPart1Forx2 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[1], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[0], equ_degree[2], 1);
				triBCdilationKernelPart1Forx2 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[1], equ_unitLength);
			}
		}
		else {
			dim3 block1(opt_degree[0], opt_degree[1], 1);
			triBCdilationKernelPart1Forx3 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[2], opt_unitLength);
			triBCdilationKernelPart1Forx3 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[2], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[0], con_degree[1], 1);
				triBCdilationKernelPart1Forx3 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[2], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[0], equ_degree[1], 1);
				triBCdilationKernelPart1Forx3 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[2], equ_unitLength);
			}
		}
	}
	else {
		if (dim == 0) {
			dim3 block1(opt_degree[1], opt_degree[2], opt_degree[3]);
			quadBCdilationKernelPart1Forx1 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[0], opt_unitLength);
			quadBCdilationKernelPart1Forx1 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[0], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[1], con_degree[2], con_degree[3]);
				quadBCdilationKernelPart1Forx1 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[0], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[1], equ_degree[2], equ_degree[3]);
				quadBCdilationKernelPart1Forx1 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[0], equ_unitLength);
			}
		}
		else if (dim == 1) {
			dim3 block1(opt_degree[0], opt_degree[2], opt_degree[3]);
			quadBCdilationKernelPart1Forx2 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[1], opt_unitLength);
			quadBCdilationKernelPart1Forx2 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[1], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[0], con_degree[2], con_degree[3]);
				quadBCdilationKernelPart1Forx2 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[1], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[0], equ_degree[2], equ_degree[3]);
				quadBCdilationKernelPart1Forx2 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[1], equ_unitLength);
			}
		}
		else if (dim == 2) {
			dim3 block1(opt_degree[0], opt_degree[1], opt_degree[3]);
			quadBCdilationKernelPart1Forx3 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[2], opt_unitLength);
			quadBCdilationKernelPart1Forx3 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[2], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[0], con_degree[1], con_degree[3]);
				quadBCdilationKernelPart1Forx3 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[2], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[0], equ_degree[1], equ_degree[3]);
				quadBCdilationKernelPart1Forx3 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[2], equ_unitLength);
			}
		}
		else {
			dim3 block1(opt_degree[0], opt_degree[1], opt_degree[2]);
			quadBCdilationKernelPart1Forx4 << < numUnit, block1 >> > (dev_opt_BC, opt_degree[3], opt_unitLength);
			quadBCdilationKernelPart1Forx4 << < numUnit * numDimension, block1 >> > (dev_pd_BC, opt_degree[3], opt_unitLength);

			if (numCons > 0) {
				dim3 block2(con_degree[0], con_degree[1], con_degree[2]);
				quadBCdilationKernelPart1Forx4 << < numUnit * numCons, block2 >> > (dev_con_BC, con_degree[3], con_unitLength);
			}

			if (numEqus > 0) {
				dim3 block3(equ_degree[0], equ_degree[1], equ_degree[2]);
				quadBCdilationKernelPart1Forx4 << < numUnit * numEqus, block3 >> > (dev_equ_BC, equ_degree[3], equ_unitLength);
			}
		}
	}

	BCdilationKernelPart2 << < numUnit, numDimension >> > (dev_interval, dev_pdFlag, dim, dev_pdValue);

	BCdilationKernelPart3 << < numUnit, numCons >> > (dev_consFlag);

	numUnit <<= 1;
}

__global__ void biBCdilationKernelPart1Forx1(float* target_BC, uint8_t x1degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x2ID = threadIdx.x;
	int BCLeftBase = unitID * unitLength + x2ID * x1degree + x1degree - 1;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + x2ID * x1degree + x1degree - 1;
	int BCRightBase = unitID * unitLength + x2ID * x1degree;
	float BCterm;
	for (int x1pos = 0; x1pos < x1degree; x1pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x1pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID] * dilationMat[x1pos][addID];
		}

		target_BC[BCLeftGrow - x1pos] = BCterm;
	}

	for (int x1pos = x1degree - 1; x1pos >= 0; x1pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x1pos; addID++) {
			BCterm += target_BC[BCRightBase + addID] * dilationMat[x1pos][addID];
		}

		target_BC[BCRightBase + x1pos] = BCterm;
	}
}

__global__ void biBCdilationKernelPart1Forx2(float* target_BC, uint8_t x2degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x1ID = threadIdx.x;
	int x1degree = blockDim.x;
	int BCLeftBase = unitID * unitLength + (x2degree - 1) * x1degree + x1ID;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + (x2degree - 1) * x1degree + x1ID;
	int BCRightBase = unitID * unitLength + x1ID;
	float BCterm;
	for (int x2pos = 0; x2pos < x2degree; x2pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x2pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID * x1degree] * dilationMat[x2pos][addID];
		}

		target_BC[BCLeftGrow - x2pos * x1degree] = BCterm;
	}

	for (int x2pos = x2degree - 1; x2pos >= 0; x2pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x2pos; addID++) {
			BCterm += target_BC[BCRightBase + addID * x1degree] * dilationMat[x2pos][addID];
		}

		target_BC[BCRightBase + x2pos * x1degree] = BCterm;
	}
}

__global__ void triBCdilationKernelPart1Forx1(float* target_BC, uint8_t x1degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x2ID = threadIdx.x;
	int x3ID = threadIdx.y;
	int x2degreePos = x1degree;
	int x3degreePos = blockDim.x * x2degreePos;
	int BCLeftBase = unitID * unitLength + x3ID * x3degreePos + x2ID * x2degreePos + x1degree - 1;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + x3ID * x3degreePos + x2ID * x2degreePos + x1degree - 1;
	int BCRightBase = unitID * unitLength + x3ID * x3degreePos + x2ID * x2degreePos;
	float BCterm;
	for (int x1pos = 0; x1pos < x1degree; x1pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x1pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID] * dilationMat[x1pos][addID];
		}

		target_BC[BCLeftGrow - x1pos] = BCterm;
	}

	for (int x1pos = x1degree - 1; x1pos >= 0; x1pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x1pos; addID++) {
			BCterm += target_BC[BCRightBase + addID] * dilationMat[x1pos][addID];
		}

		target_BC[BCRightBase + x1pos] = BCterm;
	}
}

__global__ void triBCdilationKernelPart1Forx2(float* target_BC, uint8_t x2degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x1ID = threadIdx.x;
	int x3ID = threadIdx.y;
	int x2degreePos = blockDim.x;
	int x3degreePos = x2degree * x2degreePos;
	int BCLeftBase = unitID * unitLength + x3ID * x3degreePos + (x2degree - 1) * x2degreePos + x1ID;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + x3ID * x3degreePos + (x2degree - 1) * x2degreePos + x1ID;
	int BCRightBase = unitID * unitLength + x3ID * x3degreePos + x1ID;
	float BCterm;
	for (int x2pos = 0; x2pos < x2degree; x2pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x2pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID * x2degreePos] * dilationMat[x2pos][addID];
		}

		target_BC[BCLeftGrow - x2pos * x2degreePos] = BCterm;
	}

	for (int x2pos = x2degree - 1; x2pos >= 0; x2pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x2pos; addID++) {
			BCterm += target_BC[BCRightBase + addID * x2degreePos] * dilationMat[x2pos][addID];
		}

		target_BC[BCRightBase + x2pos * x2degreePos] = BCterm;
	}
}

__global__ void triBCdilationKernelPart1Forx3(float* target_BC, uint8_t x3degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x1ID = threadIdx.x;
	int x2ID = threadIdx.y;
	int x2degreePos = blockDim.x;
	int x3degreePos = blockDim.y * x2degreePos;
	int BCLeftBase = unitID * unitLength + (x3degree - 1) * x3degreePos + x2ID * x2degreePos + x1ID;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + (x3degree - 1) * x3degreePos + x2ID * x2degreePos + x1ID;
	int BCRightBase = unitID * unitLength + x2ID * x2degreePos + x1ID;
	float BCterm;
	for (int x3pos = 0; x3pos < x3degree; x3pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x3pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID * x3degreePos] * dilationMat[x3pos][addID];
		}

		target_BC[BCLeftGrow - x3pos * x3degreePos] = BCterm;
	}

	for (int x3pos = x3degree - 1; x3pos >= 0; x3pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x3pos; addID++) {
			BCterm += target_BC[BCRightBase + addID * x3degreePos] * dilationMat[x3pos][addID];
		}

		target_BC[BCRightBase + x3pos * x3degreePos] = BCterm;
	}
}

__global__ void quadBCdilationKernelPart1Forx1(float* target_BC, uint8_t x1degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x2ID = threadIdx.x;
	int x3ID = threadIdx.y;
	int x4ID = threadIdx.z;
	int x2degreePos = x1degree;
	int x3degreePos = blockDim.x * x2degreePos;
	int x4degreePos = blockDim.y * x3degreePos;
	int BCLeftBase = unitID * unitLength + x4ID * x4degreePos + x3ID * x3degreePos + x2ID * x2degreePos + x1degree - 1;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + x4ID * x4degreePos + x3ID * x3degreePos + x2ID * x2degreePos + x1degree - 1;
	int BCRightBase = unitID * unitLength + x4ID * x4degreePos + x3ID * x3degreePos + x2ID * x2degreePos;
	float BCterm;
	for (int x1pos = 0; x1pos < x1degree; x1pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x1pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID] * dilationMat[x1pos][addID];
		}

		target_BC[BCLeftGrow - x1pos] = BCterm;
	}

	for (int x1pos = x1degree - 1; x1pos >= 0; x1pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x1pos; addID++) {
			BCterm += target_BC[BCRightBase + addID] * dilationMat[x1pos][addID];
		}

		target_BC[BCRightBase + x1pos] = BCterm;
	}
}

__global__ void quadBCdilationKernelPart1Forx2(float* target_BC, uint8_t x2degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x1ID = threadIdx.x;
	int x3ID = threadIdx.y;
	int x4ID = threadIdx.z;
	int x2degreePos = blockDim.x;
	int x3degreePos = x2degree * x2degreePos;
	int x4degreePos = blockDim.y * x3degreePos;
	int BCLeftBase = unitID * unitLength + x4ID * x4degreePos + x3ID * x3degreePos + (x2degree - 1) * x2degreePos + x1ID;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + x4ID * x4degreePos + x3ID * x3degreePos + (x2degree - 1) * x2degreePos + x1ID;
	int BCRightBase = unitID * unitLength + x4ID * x4degreePos + x3ID * x3degreePos + x1ID;
	float BCterm;
	for (int x2pos = 0; x2pos < x2degree; x2pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x2pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID * x2degreePos] * dilationMat[x2pos][addID];
		}

		target_BC[BCLeftGrow - x2pos * x2degreePos] = BCterm;
	}

	for (int x2pos = x2degree - 1; x2pos >= 0; x2pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x2pos; addID++) {
			BCterm += target_BC[BCRightBase + addID * x2degreePos] * dilationMat[x2pos][addID];
		}

		target_BC[BCRightBase + x2pos * x2degreePos] = BCterm;
	}
}

__global__ void quadBCdilationKernelPart1Forx3(float* target_BC, uint8_t x3degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x1ID = threadIdx.x;
	int x2ID = threadIdx.y;
	int x4ID = threadIdx.z;
	int x2degreePos = blockDim.x;
	int x3degreePos = blockDim.y * x2degreePos;
	int x4degreePos = x3degree * x3degreePos;
	int BCLeftBase = unitID * unitLength + x4ID * x4degreePos + (x3degree - 1) * x3degreePos + x2ID * x2degreePos + x1ID;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + x4ID * x4degreePos + (x3degree - 1) * x3degreePos + x2ID * x2degreePos + x1ID;
	int BCRightBase = unitID * unitLength + x4ID * x4degreePos + x2ID * x2degreePos + x1ID;
	float BCterm;
	for (int x3pos = 0; x3pos < x3degree; x3pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x3pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID * x3degreePos] * dilationMat[x3pos][addID];
		}

		target_BC[BCLeftGrow - x3pos * x3degreePos] = BCterm;
	}

	for (int x3pos = x3degree - 1; x3pos >= 0; x3pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x3pos; addID++) {
			BCterm += target_BC[BCRightBase + addID * x3degreePos] * dilationMat[x3pos][addID];
		}

		target_BC[BCRightBase + x3pos * x3degreePos] = BCterm;
	}
}

__global__ void quadBCdilationKernelPart1Forx4(float* target_BC, uint8_t x4degree, uint16_t unitLength)
{
	int unitID = blockIdx.x;
	int unitOffset = gridDim.x;
	int x1ID = threadIdx.x;
	int x2ID = threadIdx.y;
	int x3ID = threadIdx.z;
	int x2degreePos = blockDim.x;
	int x3degreePos = blockDim.y * x2degreePos;
	int x4degreePos = blockDim.z * x3degreePos;
	int BCLeftBase = unitID * unitLength + (x4degree - 1) * x4degreePos + x3ID * x3degreePos + x2ID * x2degreePos + x1ID;
	int BCLeftGrow = (unitID + unitOffset) * unitLength + (x4degree - 1) * x4degreePos + x3ID * x3degreePos + x2ID * x2degreePos + x1ID;
	int BCRightBase = unitID * unitLength + x3ID * x3degreePos + x2ID * x2degreePos + x1ID;
	float BCterm;
	for (int x4pos = 0; x4pos < x4degree; x4pos++) {
		BCterm = 0;

		for (int addID = 0; addID <= x4pos; addID++) {
			BCterm += target_BC[BCLeftBase - addID * x4degreePos] * dilationMat[x4pos][addID];
		}

		target_BC[BCLeftGrow - x4pos * x4degreePos] = BCterm;
	}

	for (int x4pos = x4degree - 1; x4pos >= 0; x4pos--) {
		BCterm = 0;

		for (int addID = 0; addID <= x4pos; addID++) {
			BCterm += target_BC[BCRightBase + addID * x4degreePos] * dilationMat[x4pos][addID];
		}

		target_BC[BCRightBase + x4pos * x4degreePos] = BCterm;
	}
}

__global__ void BCdilationKernelPart2(uint32_t* target_interval, bool* target_pdFlag, uint8_t dim, float* target_pdValue) {
	int unitID = blockIdx.x;
	int intID = threadIdx.x;
	int numDimension = blockDim.x;
	int numUnit = gridDim.x;

	if (intID == dim) {
		target_interval[unitID * numDimension + intID] <<= 1;
		target_interval[(unitID + numUnit) * numDimension + intID] = target_interval[unitID * numDimension + intID] + 1;
	}
	else {
		target_interval[(unitID + numUnit) * numDimension + intID] = target_interval[unitID * numDimension + intID];
	}

	target_pdValue[(unitID + numUnit) * numDimension + intID] = target_pdValue[unitID * numDimension + intID];
	target_pdFlag[(unitID + numUnit) * numDimension + intID] = target_pdFlag[unitID * numDimension + intID];
}

__global__ void BCdilationKernelPart3(char* target_consFlag) {
	int unitID = blockIdx.x;
	int valueID = threadIdx.x;
	int numCons = blockDim.x;
	int numUnit = gridDim.x;
	target_consFlag[(unitID + numUnit) * numCons + valueID] = target_consFlag[unitID * numCons + valueID];
}

void BC::findFlag() {
	if (numCons > 0) {
		BCfindFlagKernel << < numUnit * numCons, con_unitLength >> > (dev_consFlag, dev_con_BC);
		BCfindIntFlagKernel << < numUnit, 1 >> > (dev_intFlag, dev_consFlag, numCons);

		cudaMemcpy(intFlag, dev_intFlag, numUnit * sizeof(char), cudaMemcpyDeviceToHost);
	}

	if (numEqus > 0) {
		BCfindEquFlagKernel << < numUnit * numEqus, equ_unitLength >> > (dev_equsFlag, dev_equ_BC);
		BCfindEFlagKernel << < numUnit, 1 >> > (dev_eFlag, dev_equsFlag, numEqus);

		cudaMemcpy(eFlag, dev_eFlag, numUnit * sizeof(bool), cudaMemcpyDeviceToHost);
	}

	BCfindBoundKernel << < numUnit, opt_unitLength >> > (dev_bdMin, dev_bdMax, dev_opt_BC);

	cudaMemcpy(bdMin, dev_bdMin, numUnit * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(bdMax, dev_bdMax, numUnit * sizeof(float), cudaMemcpyDeviceToHost);

	BCfindDerivativeKernel << < numUnit * numDimension, opt_unitLength >> > (dev_pdFlag, dev_pd_BC, dev_pdValue);
	
	if(numDimension == 2){
		biBCfinddFlagKernel << < numUnit, 1 >> > (dev_dFlag, dev_pdFlag, numDimension, dev_interval, int_iter[0], int_iter[1]);
	}
	else if(numDimension == 3){
		triBCfinddFlagKernel << < numUnit, 1 >> > (dev_dFlag, dev_pdFlag, numDimension, dev_interval, int_iter[0], int_iter[1], int_iter[2]);
	}
	else if(numDimension == 4){
		quadBCfinddFlagKernel << < numUnit, 1 >> > (dev_dFlag, dev_pdFlag, numDimension, dev_interval, int_iter[0], int_iter[1], int_iter[2], int_iter[3]);
	}

	cudaMemcpy(dFlag, dev_dFlag, numUnit * sizeof(bool), cudaMemcpyDeviceToHost);
}

__global__ void BCfindFlagKernel(char* conFlag, float* BC)
{
	int unitID = blockIdx.x;

	if (conFlag[unitID] != 1) return;

	int tid = threadIdx.x;
	int unitLength = blockDim.x;
	__shared__ float BCbufferMin[MAX_UNIT_LENGTH];
	__shared__ float BCbufferMax[MAX_UNIT_LENGTH];

	BCbufferMin[tid] = BC[unitID * unitLength + tid];
	BCbufferMax[tid] = BC[unitID * unitLength + tid];
	__syncthreads();

	for (int i = 1; i < unitLength; i <<= 1) {
		if (tid % (i << 1) == 0 && tid + i < unitLength) {
			BCbufferMin[tid] = BCbufferMin[tid] < BCbufferMin[tid + i] ? BCbufferMin[tid] : BCbufferMin[tid + i];
			BCbufferMax[tid] = BCbufferMax[tid] > BCbufferMax[tid + i] ? BCbufferMax[tid] : BCbufferMax[tid + i];
		}
		__syncthreads();
	}

	if (tid == 0) {
		// NORMAL INEQUALITY

		if (BCbufferMin[tid] >= 0) {
			conFlag[unitID] = 0;
		}
		else if (BCbufferMax[tid] <= 0) {
			conFlag[unitID] = 2;
		}
		else {
			conFlag[unitID] = 1;
		}

		// BOUNDED INEQUALITY
		/*
		if (BCbufferMin[tid] >= 0 && BCbufferMax[tid] <= 1) {
			conFlag[unitID] = 2;
		}
		else if (BCbufferMin[tid] > 1 || BCbufferMax[tid] < 0) {
			conFlag[unitID] = 0;
		}
		else {
			conFlag[unitID] = 1;
		}
		*/
	}
}

__global__ void BCfindIntFlagKernel(char* intFlag, char* conFlag, uint8_t numCons) {
	int flagID = blockIdx.x;
	bool ifSatisfy = true;
	for (int intID = flagID * numCons; intID < (flagID + 1) * numCons; intID++) {
		if (conFlag[intID] == 0) {
			intFlag[flagID] = 0;
			return;
		}
		else if (conFlag[intID] == 1) {
			ifSatisfy = false;
		}
	}
	intFlag[flagID] = ifSatisfy ? 2 : 1;
}

__global__ void BCfindEquFlagKernel(bool* equFlag, float* BC)
{
	int unitID = blockIdx.x;

	int tid = threadIdx.x;
	int unitLength = blockDim.x;
	__shared__ float BCbufferMin[MAX_UNIT_LENGTH];
	__shared__ float BCbufferMax[MAX_UNIT_LENGTH];

	BCbufferMin[tid] = BC[unitID * unitLength + tid];
	BCbufferMax[tid] = BC[unitID * unitLength + tid];
	__syncthreads();

	for (int i = 1; i < unitLength; i <<= 1) {
		if (tid % (i << 1) == 0 && tid + i < unitLength) {
			BCbufferMin[tid] = BCbufferMin[tid] < BCbufferMin[tid + i] ? BCbufferMin[tid] : BCbufferMin[tid + i];
			BCbufferMax[tid] = BCbufferMax[tid] > BCbufferMax[tid + i] ? BCbufferMax[tid] : BCbufferMax[tid + i];
		}
		__syncthreads();
	}

	if (tid == 0) {
		if (BCbufferMin[tid] > 0 || BCbufferMax[tid] < 0) {
			equFlag[unitID] = 0;
		}
		else {
			equFlag[unitID] = 1;
		}
	}
}

__global__ void BCfindEFlagKernel(bool* eFlag, bool* equsFlag, uint8_t numEqus) {
	int flagID = blockIdx.x;
	for (int intID = flagID * numEqus; intID < (flagID + 1) * numEqus; intID++) {
		if (equsFlag[intID] == 0) {
			eFlag[flagID] = 0;
			return;
		}
	}
	eFlag[flagID] = 1;
}

__global__ void BCfindBoundKernel(float* bdMin, float* bdMax, float* BC)
{
	int unitID = blockIdx.x;
	int tid = threadIdx.x;
	int unitLength = blockDim.x;
	__shared__ float BCbufferMin[MAX_UNIT_LENGTH];
	__shared__ float BCbufferMax[MAX_UNIT_LENGTH];

	BCbufferMin[tid] = BC[unitID * unitLength + tid];
	BCbufferMax[tid] = BC[unitID * unitLength + tid];
	__syncthreads();

	for (int i = 1; i < unitLength; i <<= 1) {
		if (tid % (i << 1) == 0 && tid + i < unitLength) {
			BCbufferMin[tid] = BCbufferMin[tid] < BCbufferMin[tid + i] ? BCbufferMin[tid] : BCbufferMin[tid + i];
			BCbufferMax[tid] = BCbufferMax[tid] > BCbufferMax[tid + i] ? BCbufferMax[tid] : BCbufferMax[tid + i];
		}
		__syncthreads();
	}

	if (tid == 0) {
		bdMin[unitID] = BCbufferMin[tid];
		bdMax[unitID] = BCbufferMax[tid];
	}
}

__global__ void BCfindDerivativeKernel(bool* pdFlag, float* BC, float* pdValue)
{
	int unitID = blockIdx.x;

	if (!pdFlag[unitID]) return;

	int tid = threadIdx.x;
	int unitLength = blockDim.x;
	__shared__ float BCbufferMin[MAX_UNIT_LENGTH];
	__shared__ float BCbufferMax[MAX_UNIT_LENGTH];

	BCbufferMin[tid] = BC[unitID * unitLength + tid];
	BCbufferMax[tid] = BC[unitID * unitLength + tid];
	__syncthreads();

	for (int i = 1; i < unitLength; i <<= 1) {
		if (tid % (i << 1) == 0 && tid + i < unitLength) {
			BCbufferMin[tid] = BCbufferMin[tid] < BCbufferMin[tid + i] ? BCbufferMin[tid] : BCbufferMin[tid + i];
			BCbufferMax[tid] = BCbufferMax[tid] > BCbufferMax[tid + i] ? BCbufferMax[tid] : BCbufferMax[tid + i];
		}
		__syncthreads();
	}

	if (tid == 0) {
		if (BCbufferMin[tid] > pdValue[unitID] || BCbufferMax[tid] < pdValue[unitID]) {
			pdFlag[unitID] = false;
		}
		else {
			pdFlag[unitID] = true;
		}
	}
}

__global__ void biBCfinddFlagKernel(bool* dFlag, bool* pdFlag, uint8_t numDimension, uint32_t* interval, uint8_t iter_0, uint8_t iter_1) {
	int flagID = blockIdx.x;
	for (int pdID = flagID * numDimension; pdID < (flagID + 1) * numDimension; pdID++) {
		if (pdFlag[pdID] == false) {
			int intID = flagID * numDimension;
			if (interval[intID] == 0 || (interval[intID] + 1) == (1 << (uint32_t)iter_0)) {
				dFlag[flagID] = true;
				return;
			}
			if (interval[intID + 1] == 0 || (interval[intID + 1] + 1) == (1 << (uint32_t)iter_1)) {
				dFlag[flagID] = true;
				return;
			}

			dFlag[flagID] = false;
			return;
		}
	}
	dFlag[flagID] = true;
}

__global__ void triBCfinddFlagKernel(bool* dFlag, bool* pdFlag, uint8_t numDimension, uint32_t* interval, uint8_t iter_0, uint8_t iter_1, uint8_t iter_2) {
	int flagID = blockIdx.x;
	for (int pdID = flagID * numDimension; pdID < (flagID + 1) * numDimension; pdID++) {
		if (pdFlag[pdID] == false) {
			int intID = flagID * numDimension;
			if (interval[intID] == 0 || (interval[intID] + 1) == (1 << (uint32_t)iter_0)) {
				dFlag[flagID] = true;
				return;
			}
			if (interval[intID + 1] == 0 || (interval[intID + 1] + 1) == (1 << (uint32_t)iter_1)) {
				dFlag[flagID] = true;
				return;
			}
			if (interval[intID + 2] == 0 || (interval[intID + 2] + 1) == (1 << (uint32_t)iter_2)) {
				dFlag[flagID] = true;
				return;
			}

			dFlag[flagID] = false;
			return;
		}
	}
	dFlag[flagID] = true;
}

__global__ void quadBCfinddFlagKernel(bool* dFlag, bool* pdFlag, uint8_t numDimension, uint32_t* interval, uint8_t iter_0, uint8_t iter_1, uint8_t iter_2, uint8_t iter_3) {
	int flagID = blockIdx.x;
	for (int pdID = flagID * numDimension; pdID < (flagID + 1) * numDimension; pdID++) {
		if (pdFlag[pdID] == false) {
			int intID = flagID * numDimension;
			if (interval[intID] == 0 || (interval[intID] + 1) == (1 << (uint32_t)iter_0)) {
				dFlag[flagID] = true;
				return;
			}
			if (interval[intID + 1] == 0 || (interval[intID + 1] + 1) == (1 << (uint32_t)iter_1)) {
				dFlag[flagID] = true;
				return;
			}
			if (interval[intID + 2] == 0 || (interval[intID + 2] + 1) == (1 << (uint32_t)iter_2)) {
				dFlag[flagID] = true;
				return;
			}
			if (interval[intID + 3] == 0 || (interval[intID + 3] + 1) == (1 << (uint32_t)iter_3)) {
				dFlag[flagID] = true;
				return;
			}

			dFlag[flagID] = false;
			return;
		}
	}
	dFlag[flagID] = true;
}

void BC::eliminate() {
	estiMin = FLT_MAX;
	float estiMinMin = FLT_MAX;
	final_index = 0xffffffff;

	if (numCons > 0 && numEqus > 0) {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (eFlag[i] == 1) {
				if (intFlag[i] == 2 && bdMax[i] < estiMin) {
					estiMin = bdMax[i];
					final_index = i;
				}
				if (intFlag[i] > 0 && bdMin[i] < estiMinMin) {
					estiMinMin = bdMin[i];
				}
			}
		}
	}
	else if (numCons == 0 && numEqus > 0) {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (eFlag[i] == 1) {
				if (bdMax[i] < estiMin) {
					estiMin = bdMax[i];
					final_index = i;
				}
				if (bdMin[i] < estiMinMin) {
					estiMinMin = bdMin[i];
				}
			}
		}
	}
	else if (numCons > 0 && numEqus == 0) {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (intFlag[i] == 2 && bdMax[i] < estiMin) {
				estiMin = bdMax[i];
				final_index = i;
			}
			if (intFlag[i] > 0 && bdMin[i] < estiMinMin) {
				estiMinMin = bdMin[i];
			}
		}
	}
	else {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (bdMax[i] < estiMin) {
				estiMin = bdMax[i];
				final_index = i;
			}
			if (bdMin[i] < estiMinMin) {
				estiMinMin = bdMin[i];
			}
		}
	}

	estimated_accuracy = estiMin - estiMinMin;
	if (estimated_accuracy <= target_accuracy && estiMinMin != FLT_MAX) {
		last = true;
		return;
	}

	if ((numUnit << 1) > MAX_UNIT_NUM) return;

	uint32_t elimNum = 0;
	uint32_t saveNum = 0;
	uint32_t replaceNum;
	if (numCons > 0 && numEqus > 0) {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (intFlag[i] == 0 || eFlag[i] == 0) {
				elimPos[elimNum++] = i;
			}
			else {
				if (bdMin[i] > estiMin) {
					elimPos[elimNum++] = i;
				}
				else {
					/*
					if (dFlag[i] == false && intFlag[i] == 2) {
						elimPos[elimNum++] = i;
					}
					else {
						savePos[saveNum++] = i;
					}
					*/
					savePos[saveNum++] = i;
				}
			}
		}
	}
	else if (numCons == 0 && numEqus > 0) {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (eFlag[i] == 0) {
				elimPos[elimNum++] = i;
			}
			else {
				if (bdMin[i] > estiMin) {
					elimPos[elimNum++] = i;
				}
				else {
					savePos[saveNum++] = i;
				}
			}
		}
	}
	else if (numCons > 0 && numEqus == 0) {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (intFlag[i] == 0) {
				elimPos[elimNum++] = i;
			}
			else {
				if (bdMin[i] > estiMin) {
					elimPos[elimNum++] = i;
				}
				else {
					if (dFlag[i] == false && intFlag[i] == 2) {
						elimPos[elimNum++] = i;
					}
					else {
						savePos[saveNum++] = i;
					}
				}
			}
		}
	}
	else {
		for (uint32_t i = 0; i < numUnit; i++) {
			if (bdMin[i] > estiMin) {
				elimPos[elimNum++] = i;
			}
			else {
				if (dFlag[i] == false) {
					elimPos[elimNum++] = i;
				}
				else {
					savePos[saveNum++] = i;
				}
			}
		}
	}
	
	if (saveNum > 0 && elimNum > 0 && elimPos[0] < saveNum) {
		replaceNum = elimNum;
		for (uint32_t i = 0; i < elimNum; i++) {
			if (elimPos[i] >= saveNum) {
				replaceNum = i;
				break;
			}
		}

		cudaMemcpy(dev_elimPos, elimPos, replaceNum * sizeof(uint32_t), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_savePos, savePos, saveNum * sizeof(uint32_t), cudaMemcpyHostToDevice);

		BCeliminateKernelPart1 << < replaceNum, opt_unitLength >> > (dev_opt_BC, dev_elimPos, dev_savePos, saveNum);

		BCeliminateKernelPart2 << < replaceNum, numDimension >> > (dev_interval, dev_elimPos, dev_savePos, saveNum);

		if (numCons > 0) {
			dim3 grid1(replaceNum, numCons, 1);
			BCeliminateKernelPart3 << < grid1, con_unitLength >> > (dev_con_BC, dev_consFlag, dev_elimPos, dev_savePos, saveNum);
		}

		if (numEqus > 0) {
			BCeliminateKernelPart1 << < replaceNum * numEqus, equ_unitLength >> > (dev_equ_BC, dev_elimPos, dev_savePos, saveNum);
		}

		dim3 grid2(replaceNum, numDimension, 1);
		BCeliminateKernelPart4 << < grid2, opt_unitLength >> > (dev_pd_BC, dev_pdFlag, dev_elimPos, dev_savePos, saveNum);
	}

	numUnit = saveNum;
}

__global__ void BCeliminateKernelPart1(float* target_BC, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum) {
	int unitID = elimPos[blockIdx.x];
	int replaceID = savePos[saveNum - blockIdx.x - 1];
	int BCID = threadIdx.x;
	int unitLength = blockDim.x;
	target_BC[unitID * unitLength + BCID] = target_BC[replaceID * unitLength + BCID];
}

__global__ void BCeliminateKernelPart2(uint32_t* target_interval, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum) {
	int unitID = elimPos[blockIdx.x];
	int replaceID = savePos[saveNum - blockIdx.x - 1];
	int intID = threadIdx.x;
	int numDimension = blockDim.x;
	target_interval[unitID * numDimension + intID] = target_interval[replaceID * numDimension + intID];
}

__global__ void BCeliminateKernelPart3(float* target_BC, char* target_consFlag, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum) {
	int unitID = elimPos[blockIdx.x];
	int replaceID = savePos[saveNum - blockIdx.x - 1];
	int conID = blockIdx.y;
	int numCons = gridDim.y;
	int BCID = threadIdx.x;
	int unitLength = blockDim.x;
	target_BC[(unitID * numCons + conID) * unitLength + BCID] = target_BC[(replaceID * numCons + conID) * unitLength + BCID];

	if (BCID == 0) {
		target_consFlag[unitID * numCons + conID] = target_consFlag[replaceID * numCons + conID];
	}
}

__global__ void BCeliminateKernelPart4(float* target_BC, bool* target_pdFlag, uint32_t* elimPos, uint32_t* savePos, uint32_t saveNum) {
	int unitID = elimPos[blockIdx.x];
	int replaceID = savePos[saveNum - blockIdx.x - 1];
	int pdID = blockIdx.y;
	int numDimension = gridDim.y;
	int BCID = threadIdx.x;
	int unitLength = blockDim.x;
	target_BC[(unitID * numDimension + pdID) * unitLength + BCID] = target_BC[(replaceID * numDimension + pdID) * unitLength + BCID];

	if (BCID == 0) {
		target_pdFlag[unitID * numDimension + pdID] = target_pdFlag[replaceID * numDimension + pdID];
	}
}

void BC::finalResult() {
	final_result = new float[numDimension];
	if(final_index == 0xffffffff){
		intervalRes = new float[numUnit * numDimension];
		cudaMalloc((void**)&dev_intervalRes, numUnit * numDimension * sizeof(float));
		if (numDimension == 2) {
			biBCfinalResultKernel << < numUnit, numDimension >> > (dev_intervalRes, dev_interval, int_iter[0], int_iter[1]);
		}
		else if (numDimension == 3) {
			triBCfinalResultKernel << < numUnit, numDimension >> > (dev_intervalRes, dev_interval, int_iter[0], int_iter[1], int_iter[2]);
		}
		else {
			quadBCfinalResultKernel << < numUnit, numDimension >> > (dev_intervalRes, dev_interval, int_iter[0], int_iter[1], int_iter[2], int_iter[3]);
		}
		cudaMemcpy(intervalRes, dev_intervalRes, numUnit * numDimension * sizeof(float), cudaMemcpyDeviceToHost);

		candidates = new float[numUnit];
		estiMin = FLT_MAX;
		for (uint32_t optID = 0; optID < numUnit; optID++) {
			candidates[optID] = 0;
			for (uint16_t k = 0; k < opt->numTerms; k++)
			{
				float result = opt->coeff[k];
				for (uint16_t i = 0; i < numDimension; i++)
				{
					for (uint16_t j = 0; j < opt->degree[k * numDimension + i]; j++)
					{
						result *= intervalRes[optID * numDimension + i];
					}
				}
				candidates[optID] += result;
			}

			if (candidates[optID] < estiMin) {
				estiMin = candidates[optID];
				final_index = optID;
			}
		}

		for(int i = 0; i < numDimension; i++){
			final_result[i] = intervalRes[final_index * numDimension + i];
		}
	}
	else{
		cudaMemcpy(interval, dev_interval + final_index * numDimension, numDimension * sizeof(uint32_t), cudaMemcpyDeviceToHost);
		for(int i = 0; i < numDimension; i++){
			final_result[i] = ((float)(interval[i]) + 0.5) / (float)(1 << (uint32_t)int_iter[i]);
		}
	}
}

__global__ void biBCfinalResultKernel(float* target_intervalRes, uint32_t* interval, uint8_t iter_1, uint8_t iter_2) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (threadIdx.x == 0) {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_1));
	}
	else {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_2));
	}
}

__global__ void triBCfinalResultKernel(float* target_intervalRes, uint32_t* interval, uint8_t iter_1, uint8_t iter_2, uint8_t iter_3) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (threadIdx.x == 0) {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_1));
	}
	else if (threadIdx.x == 1) {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_2));
	}
	else {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_3));
	}
}

__global__ void quadBCfinalResultKernel(float* target_intervalRes, uint32_t* interval, uint8_t iter_1, uint8_t iter_2, uint8_t iter_3, uint8_t iter_4) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (threadIdx.x == 0) {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_1));
	}
	else if (threadIdx.x == 1) {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_2));
	}
	else if (threadIdx.x == 2) {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_3));
	}
	else {
		target_intervalRes[id] = ((float)(interval[id]) + 0.5) / (float)(1 << ((uint32_t)iter_4));
	}
}

int BC::solve(bool debugMode, bool verboseMode) {
	int exitFlag = 1; // 1 means everything is fine, -12345 means infeasible, -54321 means too many boxes have been generated before stopping criteria could be satisfied, the current result may be inaccurate
	last = false;

	findFlag();

	target_accuracy = (bdMax[0] - bdMin[0]) * STOPPING_CRITERIA;
	mexPrintf("Target accuracy: %f\n", target_accuracy);

	for (iter = 1; iter <= MAX_ITER_NUM; iter++) {
		bool ifBreak = false;
		for (dim = 0; dim < numDimension; dim++) {
			if ((numUnit << 1) > MAX_UNIT_NUM) {
				mexPrintf("Too many units, the program exits without meeting the stopping criteria!\n");
				exitFlag = -54321;
				ifBreak = true;
				break;
			}

			if (verboseMode) mexPrintf("Start iteration %d dim %d\n", iter, dim);
			dilation(dim);
			apex_numUnit = numUnit > apex_numUnit ? numUnit : apex_numUnit;
			int_iter[dim]++;
			if (verboseMode) mexPrintf("Dilation patch number: %d\n", numUnit);
			findFlag();
			if (debugMode) debug_print();
			eliminate();

			if (verboseMode) mexPrintf("Final patch number: %d\nEstimated Minimum: %.8f\nEstimated Bound: %.8f\n", numUnit, estiMin, estimated_accuracy);
			if (debugMode) debug_print();
			if (verboseMode) mexPrintf("\n");

			if (numUnit == 0) {
				mexPrintf("Infeasible!\n");
				exitFlag = -12345;
				ifBreak = true;
				break;
			}

			if (last) {
				ifBreak = true;
				break;
			}
		}
		
		mexPrintf("Finish Iteration %d\n", iter);
		mexPrintf("Final Patch number: %d\nEstimated Minimum: %.8f\n", numUnit, estiMin);

		if (ifBreak) break;
	}

	mexPrintf("Estimated accuracy: %f\n", estimated_accuracy);

	if (exitFlag != -12345) {
		finalResult();
	}

	return exitFlag;
}

#endif
