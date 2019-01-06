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
	Node* nodeList;
	Element* elementList;
	long double H, L;
	int nH, nL;
	long double k, c, ro, alpha, ambientTemperature, simulationTime, simulationStepTime, initialTemperature;
	long double** globalMatrixH;
	long double** globalMatrixHBC;
	long double** globalMatrixC;
	long double* globalVectorP;
	long double** matrixHAndVectorP;
	long double* vectorP;
	long double* vec;

	Grid();
	Grid(std::string filename);
	~Grid();

	void generateGrid(std::string filename);
	void setNodes();
	void setElements();
	void setHeatedSurfaces();
	void aggregateGrid();
	void modifyIndexes(int id, MatrixHBC tempMatrixH);
	void modifyIndexes(int id, MatrixC tempMatrixH);
	void modifyIndexes(int id, MatrixH tempMatrixH);
	void modifyIndexes(int id, VectorP tempVectorP);
	void updateVectorP();
	void calculateTemperatures();
	bool gaussMethod(int n);
	void setValueOfMatrixHAndVectorP();
	void printGridInfo();
	void printNodes();
	void printElements();
	void printGlobalMatrixH();
	void printVectorP();
};
