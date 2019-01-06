#include "Grid.h"
#include <iostream>
#include "Jacobian.h"
#include "MatrixH.h"
#include "MatrixC.h"
#include "MatrixHBC.h"
#include "VectorP.h"

int main()
{
	Grid test("test.txt");
	test.printGridInfo();
	/*
	test.printNodes();
	test.printElements();
	 */
	test.aggregateGrid();
	test.calculateTemperatures();
}
