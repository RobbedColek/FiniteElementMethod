#include "stdafx.h"
#include "Grid.h"
#include <iostream>
#include <string>

Grid::Grid()
{
}


Grid::~Grid()
{
}

void Grid::generateGrid(std::string filename)
{
	freopen(filename.c_str(), "r", stdin);

	std::cin >> H;
	std::cin >> L;
	std::cin >> nH;
	std::cin >> nL;

	std::cout << "H: " << H << std::endl;
	std::cout << "L: " << L << std::endl;
	std::cout << "nH: " << nH << std::endl;
	std::cout << "nL: " << nL << std::endl;

	// nL - liczba węzłów na długość
	// nH - liczba węzłów na wysokość

	nodeList = new Node[nH*nL];

	for (int j = 0; j < nL; j++)
	{
		for (int i = 0; i < nH; i++)
		{
			nodeList[i + (j*nH)].y = i * (H / (nH - 1));
			nodeList[i + (j*nH)].x = j * (L / (nL - 1));
		}
	}

	/*for (int i = 0; i < nH*nL; i++)
	{
		std::cout << "ID: " << i << std::endl;
		nodeList[i].printInfo();
	}*/

	elementList = new Element[(nH - 1)*(nL - 1)];

	for (int j = 0; j < (nL - 1); j++)
	{
		for (int i = 0; i < (nH - 1); i++)
		{
			elementList[i + (j*(nH - 1))].nodeID[0] = nodeList[i + (j * nH)];
			elementList[i + (j*(nH - 1))].nodeID[1] = nodeList[i + (j * nH) + nH];
			elementList[i + (j*(nH - 1))].nodeID[2] = nodeList[i + (j * nH) + nH + 1];
			elementList[i + (j*(nH - 1))].nodeID[3] = nodeList[i + (j*nH) + 1];

			elementList[i + (j*(nH - 1))].ID[0] = (i + (j * nH));
			elementList[i + (j*(nH - 1))].ID[1] = (i + (j * nH) + nH);
			elementList[i + (j*(nH - 1))].ID[2] = (i + (j * nH) + nH + 1);
			elementList[i + (j*(nH - 1))].ID[3] = (i + (j*nH) + 1);
		}
	}

	/*
	for (int i = 0; i < (nH - 1)*(nL - 1); i++)
	{
		std::cout << "ID: " << i << std::endl;
		std::cout << elementList[i].ID[0] << std::endl;
		std::cout << elementList[i].ID[1] << std::endl;
		std::cout << elementList[i].ID[2] << std::endl;
		std::cout << elementList[i].ID[3] << std::endl;
	}
	*/

	fclose(stdin);


}