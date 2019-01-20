#include "Grid.h"
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"
#include "VectorP.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <limits>
#include <iomanip>

Grid::Grid()
{
}

Grid::Grid(std::string filename)
{
	this->generateGrid(filename);
}

Grid::~Grid()
{
}

void Grid::generateGrid(std::string filename)
{
	freopen(filename.c_str(), "r", stdin);

	std::cin >> H;
	std::cin >> L;
	std::cin >> nH;
	std::cin >> nL;
	std::cin >> k;
	std::cin >> c;
	std::cin >> ro;
	std::cin >> alpha;
	std::cin >> ambientTemperature;
	std::cin >> simulationStepTime;
	std::cin >> initialTemperature;
	std::cin >> simulationTime;

	vectorP = new long double[nH];

	this->setNodes();

	this->setElements();

	this->setHeatedSurfaces();

	fclose(stdin);
}


void Grid::setNodes()
{
	nodeList = new Node[nH * nL];

	for (int j = 0; j < nL; j++)
	{
		for (int i = 0; i < nH; i++)
		{
			nodeList[i + (j * nH)].y = i * (H / (nH - 1));
			nodeList[i + (j * nH)].x = j * (L / (nL - 1));
		}
	}
}


void Grid::setElements()
{
	elementList = new Element[(nH - 1) * (nL - 1)];

	for (int j = 0; j < (nL - 1); j++)
	{
		for (int i = 0; i < (nH - 1); i++)
		{
			elementList[i + (j * (nH - 1))].nodeID[0] = nodeList[i + (j * nH)];
			elementList[i + (j * (nH - 1))].nodeID[1] = nodeList[i + (j * nH) + nH];
			elementList[i + (j * (nH - 1))].nodeID[2] = nodeList[i + (j * nH) + nH + 1];
			elementList[i + (j * (nH - 1))].nodeID[3] = nodeList[i + (j * nH) + 1];

			elementList[i + (j * (nH - 1))].ID[0] = (i + (j * nH));
			elementList[i + (j * (nH - 1))].ID[1] = (i + (j * nH) + nH);
			elementList[i + (j * (nH - 1))].ID[2] = (i + (j * nH) + nH + 1);
			elementList[i + (j * (nH - 1))].ID[3] = (i + (j * nH) + 1);
		}
	}
}

void Grid::setHeatedSurfaces()
{
	// warunki brzegowe
	for (int i = 0; i < nH - 1; i++)
	{
		elementList[i].isSurfaceHeated[3] = true;
	}

	for (int i = ((nH - 1) * (nL - 1)) - 1; i > ((nH - 1) * (nL - 1)) - 1 - (nH - 1); i--)
	{
		elementList[i].isSurfaceHeated[1] = true;
	}

	for (int i = 0; i < nL - 1; i++)
	{
		elementList[i * (nH - 1)].isSurfaceHeated[0] = true;
	}

	for (int i = 0; i < nL - 1; i++)
	{
		elementList[i * (nH - 1) + (nH - 2)].isSurfaceHeated[2] = true;
	}
}

void Grid::aggregateGrid()
{
	double temp[nL * nH];

	matrixHAndVectorP = new long double*[nL * nH];
	for (int i = 0; i < nL * nH; ++i)
		matrixHAndVectorP[i] = new long double[nL * nH];

	globalMatrixH = new long double*[nL * nH];
	for (int i = 0; i < nL * nH; ++i)
		globalMatrixH[i] = new long double[nL * nH];


	globalMatrixHBC = new long double*[nL * nH];
	for (int i = 0; i < nL * nH; ++i)
		globalMatrixHBC[i] = new long double[nL * nH];


	globalMatrixC = new long double*[nL * nH];
	for (int i = 0; i < nL * nH; ++i)
		globalMatrixC[i] = new long double[nL * nH];

	globalVectorP = new long double[nL * nH];


	for (int i = 0; i < nL * nH; i++)
	{
		globalVectorP[i] = 0;

		for (int j = 0; j < nH * nH; j++)
		{
			globalMatrixH[i][j] = 0;
			globalMatrixHBC[i][j] = 0;
			globalMatrixC[i][j] = 0;
		}
	}

	for (int i = 0; i < (nH - 1) * (nL - 1); i++)
	{
		Jacobian tempJacobian;

		tempJacobian.calculateInterpolatedCoordinates(elementList[i]);
		tempJacobian.calculateShapeFunctionsDerivatives();
		tempJacobian.calculateJacobian(elementList[i]);

		MatrixH tempMatrixH;

		tempMatrixH.calculateMatrixH(tempJacobian, k);
		modifyIndexes(i, tempMatrixH);

		MatrixC tempMatrixC;

		tempMatrixC.calculateMatrixC(tempJacobian, c, ro);
		modifyIndexes(i, tempMatrixC);

		MatrixHBC tempMatrixHBC;

		tempMatrixHBC.alpha = alpha;
		tempMatrixHBC.CalculateMatrixHBC(elementList[i]);

		modifyIndexes(i, tempMatrixHBC);

		vec = new long double[nL * nH];
		for (int j = 0; j < nL * nH; j++)
		{
			vec[j] = initialTemperature;
            temp[j] = 0;
		}

		for (int i = 0; i < nL * nH; i++)
		{
			for (int j = 0; j < nL * nH; j++)
			{
				temp[i] += ((globalMatrixC[i][j] / simulationStepTime) * vec[j]);
			}
		}

		VectorP tempVectorP;
		tempVectorP.alpha = alpha;
		tempVectorP.ambientTemperature = ambientTemperature;
		tempVectorP.CalculateVectorP(elementList[i]);

		modifyIndexes(i, tempVectorP);
	}
	for (int j = 0; j < nL * nH; j++)
	{
		vectorP[j] = globalVectorP[j];
		globalVectorP[j] = -globalVectorP[j] + temp[j];
		for (int k = 0; k < nL * nH; k++)
		{
			globalMatrixH[j][k] += globalMatrixHBC[j][k] + (globalMatrixC[j][k] / simulationStepTime);
		}
	}
}

void Grid::modifyIndexes(int id, MatrixH tempMatrixH)
{
	int tempI, tempJ;
	for (int i = 0; i < 4; i++)
	{
		switch (i)
		{
		case 0:

			tempI = elementList[id].ID[0];
			break;
		case 1:

			tempI = elementList[id].ID[1];
			break;
		case 2:

			tempI = elementList[id].ID[2];
			break;
		case 3:

			tempI = elementList[id].ID[3];
			break;
		default:

			break;
		}
		for (int j = 0; j < 4; j++)
		{
			switch (j)
			{
			case 0:

				tempJ = elementList[id].ID[0];
				break;
			case 1:

				tempJ = elementList[id].ID[1];
				break;
			case 2:

				tempJ = elementList[id].ID[2];
				break;
			case 3:

				tempJ = elementList[id].ID[3];
				break;
			default:

				break;
			}
			globalMatrixH[tempI][tempJ] += tempMatrixH.H[i][j];
		}
	}
}


void Grid::modifyIndexes(int id, MatrixHBC tempMatrixHBC)
{
	int tempI, tempJ;
	for (int i = 0; i < 4; i++)
	{
		switch (i)
		{
		case 0:

			tempI = elementList[id].ID[0];
			break;
		case 1:

			tempI = elementList[id].ID[1];
			break;
		case 2:

			tempI = elementList[id].ID[2];
			break;
		case 3:

			tempI = elementList[id].ID[3];
			break;
		default:

			break;
		}
		for (int j = 0; j < 4; j++)
		{
			switch (j)
			{
			case 0:

				tempJ = elementList[id].ID[0];
				break;
			case 1:

				tempJ = elementList[id].ID[1];
				break;
			case 2:

				tempJ = elementList[id].ID[2];
				break;
			case 3:

				tempJ = elementList[id].ID[3];
				break;
			default:

				break;
			}
			globalMatrixHBC[tempI][tempJ] += tempMatrixHBC.matrixH[i][j];
		}
	}
}


void Grid::modifyIndexes(int id, MatrixC tempMatrixC)
{
	int tempI, tempJ;
	for (int i = 0; i < 4; i++)
	{
		switch (i)
		{
		case 0:

			tempI = elementList[id].ID[0];
			break;
		case 1:

			tempI = elementList[id].ID[1];
			break;
		case 2:

			tempI = elementList[id].ID[2];
			break;
		case 3:

			tempI = elementList[id].ID[3];
			break;
		default:

			break;
		}
		for (int j = 0; j < 4; j++)
		{
			switch (j)
			{
			case 0:

				tempJ = elementList[id].ID[0];
				break;
			case 1:

				tempJ = elementList[id].ID[1];
				break;
			case 2:

				tempJ = elementList[id].ID[2];
				break;
			case 3:

				tempJ = elementList[id].ID[3];
				break;
			default:

				break;
			}
			globalMatrixC[tempI][tempJ] += tempMatrixC.C[i][j];
		}
	}
}


void Grid::modifyIndexes(int id, VectorP tempVectorP)
{
	int tempI;
	for (int i = 0; i < 4; i++)
	{
		switch (i)
		{
		case 0:

			tempI = elementList[id].ID[0];
			break;
		case 1:

			tempI = elementList[id].ID[1];
			break;
		case 2:

			tempI = elementList[id].ID[2];
			break;
		case 3:

			tempI = elementList[id].ID[3];
			break;
		default:

			break;
		}
		globalVectorP[tempI] += tempVectorP.vectorP[i];
	}
}

void Grid::setValueOfMatrixHAndVectorP()
{
	for (int i = 0; i < (nH * nL); ++i)
	{
		for (int j = 0; j < (nH * nL) + 1; ++j)
		{
			if (j == (nH * nL))
			{
				matrixHAndVectorP[i][j] = globalVectorP[i];
				continue;
			}
			matrixHAndVectorP[i][j] = globalMatrixH[i][j];
		}
	}
}

void Grid::calculateTemperatures()
{
	int n = (nH * nL);
	std::cout << "Time[s]" << std::setw(20) << "MinTemp" << std::setw(20) << "MaxTemp" << std::endl;
	for (int j = 0; j < simulationTime; j = static_cast<int>(j + simulationStepTime))
	{
        long double max = std::numeric_limits<long double>::min();
        long double min = std::numeric_limits<long double>::max();
		setValueOfMatrixHAndVectorP();
		gaussMethod(n);
		for (int i = 0; i < (nL * nH); i++)
		{
			if (vec[i] > max)
				max = vec[i];
			if (vec[i] < min)
				min = vec[i];
		}

		std::cout.setf(std::ios::fixed);
		std::cout << static_cast<int>(j + simulationStepTime) << std::setw(29) << min << std::setw(20) << max << std::endl;
		updateVectorP();
	}
}

void Grid::updateVectorP()
{
	long double val1 = 0;

	for (int j = 0; j < (nL * nH); ++j)
	{
		globalVectorP[j] = 0;
		for (int i = 0; i < (nL * nH); ++i)
		{
			val1 += ((globalMatrixC[j][i] / simulationStepTime) * vec[i]);
		}
		globalVectorP[j] = -vectorP[j] + val1;
		val1 = 0;
	}
}

// Gauss Method sourced from - https://eduinf.waw.pl/inf/alg/001_search/0076.php

bool Grid::gaussMethod(int n)
{
	const long double eps = 1e-12;
	int i, j, k;
	long double m, s;
	for (int i = 0; i < (nL * nH); i++)
	{
		vec[i] = 0;
	}

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(matrixHAndVectorP[i][i]) < eps) return false;
			m = -matrixHAndVectorP[j][i] / matrixHAndVectorP[i][i];
			for (k = i + 1; k <= n; k++)
				matrixHAndVectorP[j][k] += m * matrixHAndVectorP[i][k];
		}
	}

	for (i = n - 1; i >= 0; i--)
	{
		s = matrixHAndVectorP[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= matrixHAndVectorP[i][j] * vec[j];
		if (fabs(matrixHAndVectorP[i][i]) < eps) return false;
		vec[i] = s / matrixHAndVectorP[i][i];
	}
	return true;
}

void Grid::printGridInfo()
{
	std::cout << "---------------------------" << std::endl;
	std::cout << "Printing grid details" << std::endl;
	std::cout << "---------------------------" << std::endl;
	std::cout << "H: " << std::setw(25) << H << std::endl;
	std::cout << "L: " << std::setw(25) << L << std::endl;
	std::cout << "nH: " << std::setw(24) << nH << std::endl;
	std::cout << "nL: " << std::setw(24) << nL << std::endl;
	std::cout << "k: " << std::setw(25) << k << std::endl;
	std::cout << "c: " << std::setw(25) << c << std::endl;
	std::cout << "ro: " << std::setw(24) << ro << std::endl;
	std::cout << "alpha: " << std::setw(21) << alpha << std::endl;
	std::cout << "ambientTemperature: " << std::setw(8) << ambientTemperature << std::endl;
	std::cout << "simulationStepTime: " << std::setw(8) << simulationStepTime << std::endl;
	std::cout << "initialTemperature: " << std::setw(8) << initialTemperature << std::endl;
	std::cout << "simulationTime: " << std::setw(12) << simulationTime << std::endl;
	std::cout << "---------------------------" << std::endl;
}

void Grid::printNodes()
{
	std::cout << "---------------------------" << std::endl;
	std::cout << "Printing all nodes details" << std::endl;
	std::cout << "---------------------------" << std::endl;
	for (int i = 0; i < nH * nL; i++)
	{
		std::cout << "ID: " << i << std::endl;
		nodeList[i].printInfo();
	}
	std::cout << "---------------------------" << std::endl;
}

void Grid::printElements()
{
	std::cout << "---------------------------" << std::endl;
	std::cout << "Printing all elements details" << std::endl;
	std::cout << "---------------------------" << std::endl;
	for (int j = 0; j < (nL - 1); j++)
	{
		for (int i = 0; i < (nH - 1); i++)
		{
			elementList[i + (j * (nH - 1))].printInfo();
			elementList[i + (j * (nH - 1))].printInfo();
			elementList[i + (j * (nH - 1))].printInfo();
			elementList[i + (j * (nH - 1))].printInfo();
		}
	}
	std::cout << "---------------------------" << std::endl;
}

void Grid::printGlobalMatrixH()
{
	std::cout << "---------------------------" << std::endl;
	std::cout << "Printing Global Matrix H" << std::endl;
	std::cout << "---------------------------" << std::endl;
	for (int j = 0; j < nL * nH; j++)
	{
		for (int k = 0; k < nL * nH; k++)
		{
			std::cout << globalMatrixH[j][k] << std::endl;
		}
	}
	std::cout << "---------------------------" << std::endl;
}

void Grid::printVectorP()
{
	std::cout << "---------------------------" << std::endl;
	std::cout << "Printing Global Vector P" << std::endl;
	std::cout << "---------------------------" << std::endl;
	for (int i = 0; i < nL * nH; i++)
	{
		std::cout << globalVectorP[i] << std::endl;
	}
	std::cout << "---------------------------" << std::endl;
}
