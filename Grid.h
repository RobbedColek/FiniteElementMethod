#pragma once
#include "Node.h"
#include "Element.h"
#include <string>

class Grid
{
public:
	Node * nodeList;
	Element* elementList;
	int H, L, nH, nL;

	Grid();
	~Grid();
	void generateGrid(std::string filename);
};

