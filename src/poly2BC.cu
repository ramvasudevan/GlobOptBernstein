# ifndef POLY_2_BC_CPP_INCLUDED
# define POLY_2_BC_CPP_INCLUDED
# include "poly2BC.h"

poly::poly()
{
	numDimension = 0;
	numTerms = 0;
	coeff = NULL;
	degree = NULL;
	maxDegree = NULL;
}

poly::~poly()
{
	delete[] coeff;
	delete[] degree;
	delete[] maxDegree;
}

poly::poly(uint32_t numDimension_input, uint32_t numTerms_input, double* data_mat)
{
	numDimension = numDimension_input;
	numTerms = numTerms_input;

	coeff = new float[numTerms];
	for (uint32_t i = 0; i < numTerms; i++)
	{
		coeff[i] = (float)data_mat[i * (numDimension + 1) + numDimension];
	}

	degree = new uint32_t[numTerms * numDimension];
	maxDegree = new uint32_t[numDimension];
	for (uint32_t i = 0; i < numDimension; i++)
	{
		maxDegree[i] = 0;
		for (uint32_t j = 0; j < numTerms; j++)
		{
			degree[j * numDimension + i] = (uint32_t)data_mat[j * (numDimension + 1) + i];
			if (degree[j * numDimension + i] > maxDegree[i])
			{
				maxDegree[i] = degree[j * numDimension + i];
			}
		}
		maxDegree[i] += 1;
	}
}

poly::poly(uint32_t numDimension_input, uint32_t numTerms_input, float* data_mat)
{
	numDimension = numDimension_input;
	numTerms = numTerms_input;

	coeff = new float[numTerms];
	for (uint32_t i = 0; i < numTerms; i++)
	{
		coeff[i] = data_mat[i * (numDimension + 1) + numDimension];
	}

	degree = new uint32_t[numTerms * numDimension];
	maxDegree = new uint32_t[numDimension];
	for (uint32_t i = 0; i < numDimension; i++)
	{
		maxDegree[i] = 0;
		for (uint32_t j = 0; j < numTerms; j++)
		{
			degree[j * numDimension + i] = (uint32_t)data_mat[j * (numDimension + 1) + i];
			if (degree[j * numDimension + i] > maxDegree[i])
			{
				maxDegree[i] = degree[j * numDimension + i];
			}
		}
		maxDegree[i] += 1;
	}
}

poly::poly(uint32_t numDimension_input, uint32_t numTerms_input, double* degree_input, double* coef_input) {
	numDimension = numDimension_input;
	numTerms = numTerms_input;

	coeff = new float[numTerms];
	for (uint32_t i = 0; i < numTerms; i++)
	{
		coeff[i] = (float)coef_input[i];
	}

	degree = new uint32_t[numTerms * numDimension];
	maxDegree = new uint32_t[numDimension];
	for (uint32_t i = 0; i < numDimension; i++)
	{
		maxDegree[i] = 0;
		for (uint32_t j = 0; j < numTerms; j++)
		{
			degree[j * numDimension + i] = (uint32_t)degree_input[j * numDimension + i];
			if (degree[j * numDimension + i] > maxDegree[i])
			{
				maxDegree[i] = degree[j * numDimension + i];
			}
		}
		maxDegree[i] += 1;
	}
}

void poly::printDetails() {
	mexPrintf("DIMENSION: %d, TERMS: %d\n", numDimension, numTerms);
	for (uint32_t i = 0; i < numTerms; i++) {
		mexPrintf("%.6f,", coeff[i]);
		for (uint32_t j = 0; j < numDimension; j++) {
			mexPrintf(" %d,", degree[i * numDimension + j]);
		}
		mexPrintf("\n");
	}
}

void poly::partialDerivative(poly* &res, float &res_value, uint32_t dim) {
	float* data_mat = new float[numTerms * (numDimension + 1)];
	uint32_t pd_numTerms = 0;
	float pd_value = 0;
	for (uint32_t i = 0; i < numTerms; i++) {
		if (degree[i * numDimension + dim] > 1) {
			data_mat[pd_numTerms * (numDimension + 1) + numDimension] = coeff[i] * degree[i * numDimension + dim];

			for (uint32_t j = 0; j < numDimension; j++) {
				data_mat[pd_numTerms * (numDimension + 1) + j] = degree[i * numDimension + j];
			}
			data_mat[pd_numTerms * (numDimension + 1) + dim]--;

			pd_numTerms++;
		}
		else if (degree[i * numDimension + dim] == 1) {
			uint32_t temp = 0;
			for (uint32_t j = 0; j < numDimension; j++) {
				if (j != dim) {
					temp += degree[i * numDimension + j];
				}
			}

			if (temp > 0) {
				data_mat[pd_numTerms * (numDimension + 1) + numDimension] = coeff[i] * degree[i * numDimension + dim];

				for (uint32_t j = 0; j < numDimension; j++) {
					data_mat[pd_numTerms * (numDimension + 1) + j] = degree[i * numDimension + j];
				}
				data_mat[pd_numTerms * (numDimension + 1) + dim]--;

				pd_numTerms++;
			}
			else {
				pd_value -= coeff[i];
			}
		}
	}

	res = new poly(numDimension, pd_numTerms, data_mat);
	res_value = pd_value;

	delete[] data_mat;
}

# endif