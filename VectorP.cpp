//
// Created by samsung on 1/3/19.
//
#include "VectorP.h"
#include <cmath>
#include <iostream>

VectorP::VectorP()
{
	Px[0] = -1 / sqrt(3);
	Py[0] = -1;

	Px[1] = 1 / sqrt(3);
	Py[1] = -1;

	Px[2] = 1;
	Py[2] = -1 / sqrt(3);

	Px[3] = 1;
	Py[3] = 1 / sqrt(3);

	Px[4] = 1 / sqrt(3);
	Py[4] = 1;

	Px[5] = -1 / sqrt(3);
	Py[5] = 1;

	Px[6] = -1;
	Py[6] = 1 / sqrt(3);

	Px[7] = -1;
	Py[7] = -1 / sqrt(3);
}

VectorP::~VectorP()
{
}

void VectorP::CalculateVectorP(Element element)
{
	long double sum1[4], sum2[4], sum3[4], sum4[4];

	long double length[4], detJ[4];

	length[0] = sqrt(
		pow(element.nodeID[1].x - element.nodeID[0].x, 2) + pow(element.nodeID[1].y - element.nodeID[0].y, 2));
	length[1] = sqrt(
		pow(element.nodeID[1].x - element.nodeID[2].x, 2) + pow(element.nodeID[1].y - element.nodeID[2].y, 2));
	length[2] = sqrt(
		pow(element.nodeID[2].x - element.nodeID[3].x, 2) + pow(element.nodeID[2].y - element.nodeID[3].y, 2));
	length[3] = sqrt(
		pow(element.nodeID[0].x - element.nodeID[3].x, 2) + pow(element.nodeID[0].y - element.nodeID[3].y, 2));

	detJ[0] = length[0] / 2;
	detJ[1] = length[1] / 2;
	detJ[2] = length[2] / 2;
	detJ[3] = length[3] / 2;

	sum1[0] = ((0.25 * (1 - Px[0]) * (1 - Py[0])) + (0.25 * (1 - Px[1]) * (1 - Py[1]))) * detJ[0];
	sum1[1] = ((0.25 * (1 + Px[0]) * (1 - Py[0])) + (0.25 * (1 + Px[1]) * (1 - Py[1]))) * detJ[0];
	sum1[2] = ((0.25 * (1 + Px[0]) * (1 + Py[0])) + (0.25 * (1 + Px[1]) * (1 + Py[1]))) * detJ[0];
	sum1[3] = ((0.25 * (1 - Px[0]) * (1 + Py[0])) + (0.25 * (1 - Px[1]) * (1 + Py[1]))) * detJ[0];

	sum2[0] = ((0.25 * (1 - Px[2]) * (1 - Py[2])) + (0.25 * (1 - Px[3]) * (1 - Py[3]))) * detJ[1];
	sum2[1] = ((0.25 * (1 + Px[2]) * (1 - Py[2])) + (0.25 * (1 + Px[3]) * (1 - Py[3]))) * detJ[1];
	sum2[2] = ((0.25 * (1 + Px[2]) * (1 + Py[2])) + (0.25 * (1 + Px[3]) * (1 + Py[3]))) * detJ[1];
	sum2[3] = ((0.25 * (1 - Px[2]) * (1 + Py[2])) + (0.25 * (1 - Px[3]) * (1 + Py[3]))) * detJ[1];

	sum3[0] = ((0.25 * (1 - Px[4]) * (1 - Py[4])) + (0.25 * (1 - Px[5]) * (1 - Py[5]))) * detJ[2];
	sum3[1] = ((0.25 * (1 + Px[4]) * (1 - Py[4])) + (0.25 * (1 + Px[5]) * (1 - Py[5]))) * detJ[2];
	sum3[2] = ((0.25 * (1 + Px[4]) * (1 + Py[4])) + (0.25 * (1 + Px[5]) * (1 + Py[5]))) * detJ[2];
	sum3[3] = ((0.25 * (1 - Px[4]) * (1 + Py[4])) + (0.25 * (1 - Px[5]) * (1 + Py[5]))) * detJ[2];

	sum4[0] = ((0.25 * (1 - Px[6]) * (1 - Py[6])) + (0.25 * (1 - Px[7]) * (1 - Py[7]))) * detJ[3];
	sum4[1] = ((0.25 * (1 + Px[6]) * (1 - Py[6])) + (0.25 * (1 + Px[7]) * (1 - Py[7]))) * detJ[3];
	sum4[2] = ((0.25 * (1 + Px[6]) * (1 + Py[6])) + (0.25 * (1 + Px[7]) * (1 + Py[7]))) * detJ[3];
	sum4[3] = ((0.25 * (1 - Px[6]) * (1 + Py[6])) + (0.25 * (1 - Px[7]) * (1 + Py[7]))) * detJ[3];


	for (int i = 0; i < 4; i++)
	{
		vectorP[i] = (-1) * ((element.isSurfaceHeated[0] * sum1[i] + element.isSurfaceHeated[1] * sum2[i] + element.
			isSurfaceHeated[2] * sum3[i] + element.isSurfaceHeated[3] * sum4[i]) * ambientTemperature * alpha);
	}
}
