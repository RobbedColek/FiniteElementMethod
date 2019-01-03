#pragma once
#include "Node.h"
#include "Element.h"
#include "Jacobian.h"
#include <string>

class Grid
{
public:
	Node * nodeList;
	Element* elementList;
	double **globalMatrixH;
	float H, L;
	int nH, nL;
	int tempI, tempJ;

	Grid();
	~Grid();
	void generateGrid(std::string filename);
	void aggregateGrid(double k, double c, double ro, double alpha, double ambientTemperature);
};

