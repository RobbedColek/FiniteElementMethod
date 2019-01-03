//
// Created by samsung on 1/3/19.
//

#ifndef FINITEELEMENTMETHOD_VECTORP_H
#define FINITEELEMENTMETHOD_VECTORP_H

#include "Element.h"


class VectorP {
public:
    double alpha;
    double Px[8], Py[8];
    double vectorP[4];
    double ambientTemperature;

    VectorP();
    ~VectorP();

    void CalculateVectorP(Element element);
};


#endif //FINITEELEMENTMETHOD_VECTORP_H
