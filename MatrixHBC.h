//
// Created by samsung on 1/3/19.
//

#ifndef FINITEELEMENTMETHOD_MATRIXHBC_H
#define FINITEELEMENTMETHOD_MATRIXHBC_H


#include "Element.h"

class MatrixHBC {
public:
    long double alpha;
    long double Px[8], Py[8];
    long double matrixH[4][4];

    MatrixHBC();
    ~MatrixHBC();

    void CalculateMatrixHBC(Element element);

};


#endif //FINITEELEMENTMETHOD_MATRIXHBC_H
