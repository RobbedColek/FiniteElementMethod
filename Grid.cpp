#include "Grid.h"
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"
#include "VectorP.h"
#include <iostream>
#include <string>

Grid::Grid()
{
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

	std::cout << "H: " << H << std::endl;
	std::cout << "L: " << L << std::endl;
	std::cout << "nH: " << nH << std::endl;
	std::cout << "nL: " << nL << std::endl;

	// nL - liczba w?z??w na d?ugo??
	// nH - liczba w?z??w na wysoko??

	nodeList = new Node[nH*nL];

	for (int j = 0; j < nL; j++)
	{
		for (int i = 0; i < nH; i++)
		{
			nodeList[i + (j*nH)].y = i * (H / (nH - 1));
			nodeList[i + (j*nH)].x = j * (L / (nL - 1));
		}
	}

	/*for (int i = 0; i < nH*nL; i++)
	{
		std::cout << "ID: " << i << std::endl;
		nodeList[i].printInfo();
	}*/

	elementList = new Element[(nH - 1)*(nL - 1)];

	for (int j = 0; j < (nL - 1); j++)
	{
		for (int i = 0; i < (nH - 1); i++)
		{
			elementList[i + (j*(nH - 1))].nodeID[0] = nodeList[i + (j * nH)];
			elementList[i + (j*(nH - 1))].nodeID[1] = nodeList[i + (j * nH) + nH];
			elementList[i + (j*(nH - 1))].nodeID[2] = nodeList[i + (j * nH) + nH + 1];
			elementList[i + (j*(nH - 1))].nodeID[3] = nodeList[i + (j*nH) + 1];

			elementList[i + (j*(nH - 1))].ID[0] = (i + (j * nH));
			elementList[i + (j*(nH - 1))].ID[1] = (i + (j * nH) + nH);
			elementList[i + (j*(nH - 1))].ID[2] = (i + (j * nH) + nH + 1);
			elementList[i + (j*(nH - 1))].ID[3] = (i + (j*nH) + 1);
		}
	}

	// warunki brzegowe
	for(int i = 0; i < nH-1; i++)
	{
		elementList[i].isSurfaceHeated[3] = true;
	}

	for(int i = ((nH-1)*(nL-1))-1; i > ((nH-1)*(nL-1))-1-(nH-1); i--)
	{
		elementList[i].isSurfaceHeated[1] = true;
	}

	for(int i = 0; i < nL-1; i++)
	{
		elementList[i*(nH-1)].isSurfaceHeated[0] = true;
	}

	for(int i = 0; i < nL-1; i++)
	{
		elementList[i*(nH-1)+(nH-2)].isSurfaceHeated[2] = true;
	}

	/*
	for (int i = 0; i < (nH - 1)*(nL - 1); i++)
	{
	    std::cout << "ID: " << i << std::endl;
	    std::cout << elementList[i].isSurfaceHeated[0] << std::endl;
	    std::cout << elementList[i].isSurfaceHeated[1] << std::endl;
	    std::cout << elementList[i].isSurfaceHeated[2] << std::endl;
	    std::cout << elementList[i].isSurfaceHeated[3] << std::endl;
	}
	*/

	fclose(stdin);


}

void Grid::aggregateGrid(double k, double c, double ro, double alpha, double ambientTemperature) {
	// k - conductivity (przewodnosc)
	// c - cieplo wlasciwe
	// ro - gestosc
	// alpha - convection (konwekcja)
	// ambientTemperature - temp otoczenia

	double** globalMatrixH = new double*[nL*nH];
	for(int i = 0; i < nL*nH; ++i)
		globalMatrixH[i] = new double[nL*nH];

	for(int i = 0; i < nL*nH; i++)
	{
		for(int j = 0; j < nH*nH; j++)
		{
			globalMatrixH[i][j] = 0;
		}
	}

	for (int i = 0; i < (nL-1)*(nH-1); i++)
	{
		for (int z = 0; z < 4; z++)
		{
			for (int l = 0; l <4 ; l++)
			{
				int indexOuter = elementList[i].ID[z]-1;
				int indexInner = elementList[i].ID[l]-1;

				Jacobian tempJacobian;

				tempJacobian.calculateInterpolatedCoordinates(elementList[i]);
				tempJacobian.calculateShapeFunctionsDerivatives();
				tempJacobian.calculateJacobian(elementList[i]);

				MatrixH tempMatrixH;

				tempMatrixH.calculateMatrixH(tempJacobian, k);

				globalMatrixH[indexOuter][indexInner] = tempMatrixH.H[l][z];

				std::cout << globalMatrixH[indexOuter][indexInner] << std::endl;

				//globalMatrixC[indexOuter][indexInner] += element[i].tableMatrixC[l][k];
				//globalVectorP[indexOuter][indexInner] += element[i].tableVectorP[l][k];
			}
		}
	}

	/*
	for(int i = 0; i < 1; i++) {
			Jacobian tempJacobian;

			tempJacobian.calculateInterpolatedCoordinates(elementList[i]);
			tempJacobian.calculateShapeFunctionsDerivatives();
			tempJacobian.calculateJacobian(elementList[i]);

			MatrixH tempMatrixH;

			tempMatrixH.calculateMatrixH(tempJacobian, k);

			MatrixC tempMatrixC;

			tempMatrixC.calculateMatrixC(tempJacobian, c, ro);

			MatrixHBC tempMatrixHBC;

			tempMatrixHBC.alpha = alpha;

			tempMatrixHBC.CalculateMatrixHBC(elementList[i]);

			VectorP tempVectorP;
			tempVectorP.alpha = alpha;
			tempVectorP.ambientTemperature = ambientTemperature;
			tempVectorP.CalculateVectorP(elementList[i]);


			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					matrixH[j][k] = matrixH[j][k] + tempMatrixH.H[j][k];
					std::cout << matrixH[j][k] << std::endl;
				}
			}
		}
	 */
}
