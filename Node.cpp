#include "stdafx.h"
#include "Node.h"
#include <iostream>


Node::Node()
{

}

Node::Node(double xx, double yy, double tt)
{
	x = xx;
	y = yy;
	t = tt;
}

Node::~Node()
{
}

void Node::printInfo()
{
	std::cout << x << std::endl;
	std::cout << y << std::endl;
}
