#include "stdafx.h"
#include "MatrixH.h"
#include <iostream>

MatrixH::MatrixH()
{
}

MatrixH::~MatrixH()
{
}

void MatrixH::calculateMatrixH(Jacobian jacobian, double k)
{
	for(int i = 0; i < 4; i++)
	{
		dNdx[0][i] = jacobian.inversedJacobian[0][0] * jacobian.dNdksi[i][0] + jacobian.inversedJacobian[0][1] * jacobian.dNdeta[i][0];
		dNdx[1][i] = jacobian.inversedJacobian[1][0] * jacobian.dNdksi[i][1] + jacobian.inversedJacobian[1][1] * jacobian.dNdeta[i][1];
		dNdx[2][i] = jacobian.inversedJacobian[2][0] * jacobian.dNdksi[i][2] + jacobian.inversedJacobian[2][1] * jacobian.dNdeta[i][2];
		dNdx[3][i] = jacobian.inversedJacobian[3][0] * jacobian.dNdksi[i][3] + jacobian.inversedJacobian[3][1] * jacobian.dNdeta[i][3];

		dNdy[0][i] = jacobian.inversedJacobian[0][2] * jacobian.dNdksi[i][0] + jacobian.inversedJacobian[0][3] * jacobian.dNdeta[i][0];
		dNdy[1][i] = jacobian.inversedJacobian[1][2] * jacobian.dNdksi[i][1] + jacobian.inversedJacobian[1][3] * jacobian.dNdeta[i][1];
		dNdy[2][i] = jacobian.inversedJacobian[2][2] * jacobian.dNdksi[i][2] + jacobian.inversedJacobian[2][3] * jacobian.dNdeta[i][2];
		dNdy[3][i] = jacobian.inversedJacobian[3][2] * jacobian.dNdksi[i][3] + jacobian.inversedJacobian[3][3] * jacobian.dNdeta[i][3];
	}
	//[1][2][3] 1-nr tabeli 2-nr wiersza 3-nr kolumny
	
	for(int i=0;i<4;i++)
	{
		for (int j = 0; j<4; j++)
		{
			for (int z = 0; z<4; z++)
			{
				dNdxT[i][z][j] = dNdx[i][j] * dNdx[i][z];
				dNdyT[i][z][j] = dNdy[i][j] * dNdy[i][z];

				dNdxTdetJ[i][z][j] = dNdxT[i][z][j]*jacobian.detJacobian[i];
				dNdyTdetJ[i][z][j] = dNdyT[i][z][j]*jacobian.detJacobian[i];

				sum[i][z][j] = (dNdxTdetJ[i][z][j] + dNdyTdetJ[i][z][j])*k;
			}
		}
	}
	for (int i = 0; i<4; i++)
	{
		for (int j = 0; j<4; j++)
		{
			H[j][i] = sum[0][j][i] + sum[1][j][i] + sum[2][j][i] + sum[3][j][i];
			std::cout << H[j][i] << std::endl;
		}
	}
}