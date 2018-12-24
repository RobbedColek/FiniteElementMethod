#include "Grid.h"
#include <iostream>
#include "Jacobian.h"
#include "MatrixH.h"

int main()
{
	Grid test;
	test.generateGrid("test.txt");

	Jacobian testJacobian;

	Element testElement;

	testElement.nodeID[0].x = 0;
	testElement.nodeID[0].y = 0;

	testElement.nodeID[1].x = 0.025;
	testElement.nodeID[1].y = 0;

	testElement.nodeID[2].x = 0.025;
	testElement.nodeID[2].y = 0.025;

	testElement.nodeID[3].x = 0;
	testElement.nodeID[3].y = 0.025;

	testJacobian.calculateInterpolatedCoordinates(testElement);

	testJacobian.calculateShapeFunctionsDerivatives();

	testJacobian.calculateJacobian(testElement);

	MatrixH testMatrixH;

	testMatrixH.calculateMatrixH(testJacobian, 30);
}
