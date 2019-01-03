#pragma once
#include "Node.h"
#include "Element.h"
#include "MatrixH.h"
#include "MatrixHBC.h"
#include "MatrixC.h"
#include "VectorP.h"
#include <string>

class Grid
{
public:
	Node * nodeList;
	Element* elementList;
	double **globalMatrixH;
	double **globalMatrixHBC;
	double **globalMatrixC;
	double *globalVectorP;
	float H, L;
	int nH, nL;
	double k, c, ro, alpha, ambientTemperature, simulationTime, simulationStepTime, initialTemperature;
	double *vec;
	double **matrixHAndVectorP;
	double *vectorP;

	Grid();
	~Grid();
	void generateGrid(std::string filename);
	void aggregateGrid();
	void modifyIndexes(int id, MatrixHBC tempMatrixH);
	void modifyIndexes(int id, MatrixC tempMatrixH);
	void modifyIndexes(int id, MatrixH tempMatrixH);
	void modifyIndexes(int id, VectorP tempVectorP);
	void updateVectorP();
	void calculateTemperatures();
	bool gaussMethod(int n);
	void setValueOfMatrixHAndVectorP();
};
