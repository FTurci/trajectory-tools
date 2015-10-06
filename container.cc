#include "container.h"
#include <iostream>
Container::Container()
{
}

Container::Container(const Container& copy) : boundaries(copy.boundaries)
{
    
}