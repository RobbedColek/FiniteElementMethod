//
// Created by samsung on 1/3/19.
//

#ifndef FINITEELEMENTMETHOD_VECTORP_H
#define FINITEELEMENTMETHOD_VECTORP_H

#include "Element.h"

class VectorP
{
public:
	long double alpha;
	long double Px[8], Py[8];
	long double vectorP[4];
	long double ambientTemperature;

	VectorP();
	~VectorP();

	void CalculateVectorP(Element element);
};


#endif //FINITEELEMENTMETHOD_VECTORP_H
