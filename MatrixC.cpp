#include "MatrixC.h"
#include <iostream>

MatrixC::MatrixC()
{
}

MatrixC::~MatrixC()
{
}

void MatrixC::calculateMatrixC(Jacobian jacobian, double c, double ro)
{
    double tempMatrix[4][4][4];

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            tempMatrix[i][j][0] = jacobian.N[0][j] * jacobian.N[0][i] * jacobian.detJacobian[0] * c * ro;
            tempMatrix[i][j][1] = jacobian.N[1][j] * jacobian.N[1][i] * jacobian.detJacobian[1] * c * ro;
            tempMatrix[i][j][2] = jacobian.N[2][j] * jacobian.N[2][i] * jacobian.detJacobian[2] * c * ro;
            tempMatrix[i][j][3] = jacobian.N[3][j] * jacobian.N[3][i] * jacobian.detJacobian[3] * c * ro;

            C[i][j] = tempMatrix[i][j][0] + tempMatrix[i][j][1] + tempMatrix[i][j][2] + tempMatrix[i][j][3];

            // std::cout << C[i][j] << std::endl;
        }
    }
}