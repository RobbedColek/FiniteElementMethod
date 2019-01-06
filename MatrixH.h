#pragma once
#include "Jacobian.h"

class MatrixH
{
public:
	long double H[4][4], dNdx[4][4], dNdy[4][4], dNdxT[4][4][4], dNdyT[4][4][4], dNdxTdetJ[4][4][4], dNdyTdetJ[4][4][4],
	            sum[4][4][4];

	MatrixH();
	~MatrixH();

	void calculateMatrixH(Jacobian jacobian, long double k);
};
