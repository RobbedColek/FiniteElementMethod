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
