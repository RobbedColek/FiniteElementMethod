#pragma once
#include "Jacobian.h"

class MatrixC
{
public:
	long double C[4][4], dNdx[4][4], dNdy[4][4], dNdxT[4][4][4], dNdyT[4][4][4], dNdxTdetJ[4][4][4], dNdyTdetJ[4][4][4],
	            sum[4][4][4];

	MatrixC();
	~MatrixC();

	void calculateMatrixC(Jacobian jacobian, long double c, long double ro);
};
