#include <iostream>
#include "Element.h"



Element::Element()
{
    for(int i = 0; i < 4; i++)
    {
        isSurfaceHeated[i] = false;
    }
}


Element::~Element()
{
}

void Element::printInfo() {
    for(int i = 0; i < 4; i++)
    {
        std::cout << "Node ID in this element: " << ID[i] << std::endl;
    }
    std::cout << "Is heated from bottom: " << isSurfaceHeated[0] << std::endl;
    std::cout << "Is heated from right: " << isSurfaceHeated[1] << std::endl;
    std::cout << "Is heated from top: " << isSurfaceHeated[2] << std::endl;
    std::cout << "Is heated from left: " << isSurfaceHeated[3] << std::endl;
}
