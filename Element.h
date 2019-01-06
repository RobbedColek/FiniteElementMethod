#pragma once
#include "Node.h"

class Element
{
public:
	Node nodeID[4];  
	int ID[4];
	bool isSurfaceHeated[4];

	Element();
	~Element();
	void printInfo();
};