//
// Created by samsung on 1/3/19.
//

#include <cmath>
#include <iostream>
#include "MatrixHBC.h"

MatrixHBC::MatrixHBC() {
    Px[0] = -1/sqrt(3);
    Py[0] = -1;

    Px[1] = 1/sqrt(3);
    Py[1] = -1;

    Px[2] = 1;
    Py[2] = -1/sqrt(3);

    Px[3] = 1;
    Py[3] = 1/sqrt(3);

    Px[4] = 1/sqrt(3);
    Py[4] = 1;

    Px[5] = -1/sqrt(3);
    Py[5] = 1;

    Px[6] = -1;
    Py[6] = 1/sqrt(3);

    Px[7] = -1;
    Py[7] = -1/sqrt(3);
}

MatrixHBC::~MatrixHBC() {

}

void MatrixHBC::CalculateMatrixHBC(Element element) {
    long double sum1[4][4], sum2[4][4], sum3[4][4], sum4[4][4];

    long double temp1[4], temp2[4];

    long double length[4], detJ[4];

    length[0] = sqrt(pow(element.nodeID[1].x-element.nodeID[0].x, 2) + pow(element.nodeID[1].y - element.nodeID[0].y, 2));
    length[1] = sqrt(pow(element.nodeID[1].x-element.nodeID[2].x, 2) + pow(element.nodeID[1].y - element.nodeID[2].y, 2));
    length[2] = sqrt(pow(element.nodeID[2].x-element.nodeID[3].x, 2) + pow(element.nodeID[2].y - element.nodeID[3].y, 2));
    length[3] = sqrt(pow(element.nodeID[0].x-element.nodeID[3].x, 2) + pow(element.nodeID[0].y - element.nodeID[3].y, 2));

    detJ[0] = length[0]/2;
    detJ[1] = length[1]/2;
    detJ[2] = length[2]/2;
    detJ[3] = length[3]/2;

    temp1[0] = (0.25 * (1 - Px[0]) * (1 - Py[0]));
    temp1[1] = (0.25 * (1 + Px[0]) * (1 - Py[0]));
    temp1[2] = (0.25 * (1 + Px[0]) * (1 + Py[0]));
    temp1[3] = (0.25 * (1 - Px[0]) * (1 + Py[0]));

    temp2[0] = (0.25 * (1 - Px[1]) * (1 - Py[1]));
    temp2[1] = (0.25 * (1 + Px[1]) * (1 - Py[1]));
    temp2[2] = (0.25 * (1 + Px[1]) * (1 + Py[1]));
    temp2[3] = (0.25 * (1 - Px[1]) * (1 + Py[1]));


        for(int i = 0; i < 4; i++)
        {
                for(int j = 0; j < 4; j++)
                {
                        sum1[i][j] = ((alpha * temp1[i] * temp1[j]) + (alpha * temp2[i] * temp2[j])) * detJ[0];

                        //std::cout << sum1[i][j] << std::endl;
                }
        }

        temp1[0] = (0.25 * (1 - Px[2]) * (1 - Py[2]));
        temp1[1] = (0.25 * (1 + Px[2]) * (1 - Py[2]));
        temp1[2] = (0.25 * (1 + Px[2]) * (1 + Py[2]));
        temp1[3] = (0.25 * (1 - Px[2]) * (1 + Py[2]));

        temp2[0] = (0.25 * (1 - Px[3]) * (1 - Py[3]));
        temp2[1] = (0.25 * (1 + Px[3]) * (1 - Py[3]));
        temp2[2] = (0.25 * (1 + Px[3]) * (1 + Py[3]));
        temp2[3] = (0.25 * (1 - Px[3]) * (1 + Py[3]));


        for(int i = 0; i < 4; i++)
        {
                for(int j = 0; j < 4; j++)
                {
                        sum2[i][j] = ((alpha * temp1[i] * temp1[j]) + (alpha * temp2[i] * temp2[j])) * detJ[1];

                        //std::cout << sum2[i][j] << std::endl;
                }
        }

        temp1[0] = (0.25 * (1 - Px[4]) * (1 - Py[4]));
        temp1[1] = (0.25 * (1 + Px[4]) * (1 - Py[4]));
        temp1[2] = (0.25 * (1 + Px[4]) * (1 + Py[4]));
        temp1[3] = (0.25 * (1 - Px[4]) * (1 + Py[4]));

        temp2[0] = (0.25 * (1 - Px[5]) * (1 - Py[5]));
        temp2[1] = (0.25 * (1 + Px[5]) * (1 - Py[5]));
        temp2[2] = (0.25 * (1 + Px[5]) * (1 + Py[5]));
        temp2[3] = (0.25 * (1 - Px[5]) * (1 + Py[5]));


        for(int i = 0; i < 4; i++)
        {
                for(int j = 0; j < 4; j++)
                {
                        sum3[i][j] = ((alpha * temp1[i] * temp1[j]) + (alpha * temp2[i] * temp2[j])) * detJ[2];

                        //std::cout << sum3[i][j] << std::endl;
                }
        }

        temp1[0] = (0.25 * (1 - Px[6]) * (1 - Py[6]));
        temp1[1] = (0.25 * (1 + Px[6]) * (1 - Py[6]));
        temp1[2] = (0.25 * (1 + Px[6]) * (1 + Py[6]));
        temp1[3] = (0.25 * (1 - Px[6]) * (1 + Py[6]));

        temp2[0] = (0.25 * (1 - Px[7]) * (1 - Py[7]));
        temp2[1] = (0.25 * (1 + Px[7]) * (1 - Py[7]));
        temp2[2] = (0.25 * (1 + Px[7]) * (1 + Py[7]));
        temp2[3] = (0.25 * (1 - Px[7]) * (1 + Py[7]));


        for(int i = 0; i < 4; i++)
        {
                for(int j = 0; j < 4; j++)
                {
                        sum4[i][j] = ((alpha * temp1[i] * temp1[j]) + (alpha * temp2[i] * temp2[j])) * detJ[3];

                        //std::cout << sum4[i][j] << std::endl;
                }
        }

        for(int i = 0; i < 4; i++)
        {
                for(int j = 0; j < 4; j++)
                {
                        matrixH[i][j] = element.isSurfaceHeated[0] * sum1[i][j] + element.isSurfaceHeated[1] * sum2[i][j] + element.isSurfaceHeated[2] * sum3[i][j] + element.isSurfaceHeated[3] * sum4[i][j];

                        // std::cout << matrixH[i][j] << std::endl;
                }
        }

}

