#ifndef __CONTAINER_H
#define __CONTAINER_H

#include <vector>


class Container
{
public:
    Container();

protected:
    inline double apply_boundaries(double value, unsigned int dimension)
    {
        if(value>dimension*0.5) return value-dimension;
        else if(value<-dimension*0.5) return value+dimension;
        else return value;
    }
    std::vector<double> boundaries;
};

#endif
