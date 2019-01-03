#include "Grid.h"
#include <iostream>
#include "Jacobian.h"
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"
#include "VectorP.h"

int main()
{
	Grid test;
	test.generateGrid("test.txt");

	/*
	Jacobian testJacobian;

	Element testElement;

	testElement.nodeID[0].x = 0;
	testElement.nodeID[0].y = 0;

	testElement.nodeID[1].x = 0.03333333333333333333333;
	testElement.nodeID[1].y = 0;

	testElement.nodeID[2].x = 0.03333333333333333333333;
	testElement.nodeID[2].y = 0.03333333333333333333333;

	testElement.nodeID[3].x = 0;
	testElement.nodeID[3].y = 0.03333333333333333333333;

	testElement.isSurfaceHeated[0] = true;
	testElement.isSurfaceHeated[3] = true;

	testJacobian.calculateInterpolatedCoordinates(testElement);

	testJacobian.calculateShapeFunctionsDerivatives();

	testJacobian.calculateJacobian(testElement);

	MatrixH testMatrixH;

	testMatrixH.calculateMatrixH(testJacobian, 25);

	MatrixC testMatrixC;

	testMatrixC.calculateMatrixC(testJacobian, 700, 7800);

	MatrixHBC testMatrixHBC;

	testMatrixHBC.alpha = 300;

	testMatrixHBC.CalculateMatrixHBC(testElement);

	 VectorP testVectorP;
	 testVectorP.alpha = 300;
	 testVectorP.ambientTemperature = 1200;
	 testVectorP.CalculateVectorP(testElement);
	 */

	test.aggregateGrid(25, 700, 7800, 300, 1200);
}
