#pragma once
#include "Element.h"

class MatrixHBC
{
public:
	long double alpha;
	long double Px[8], Py[8];
	long double matrixH[4][4];

	MatrixHBC();
	~MatrixHBC();

	void CalculateMatrixHBC(Element element);
};