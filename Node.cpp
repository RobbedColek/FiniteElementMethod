#include "Node.h"
#include <iostream>


Node::Node()
{

}

Node::~Node()
{
}

void Node::printInfo()
{
	std::cout << "X: " << x << std::endl;
	std::cout << "Y: " << y << std::endl;
}
