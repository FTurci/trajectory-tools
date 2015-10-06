#ifndef __CONTAINER_H
#define __CONTAINER_H

#include <vector>


class Container
{
public:
    Container();
    Container(const Container& copy);

protected:
    inline double apply_boundaries(double value, unsigned int dimension)
    {
        if(value > this->boundaries[dimension]*0.5) return value-this->boundaries[dimension];
        else if(value < -this->boundaries[dimension]*0.5) return value+this->boundaries[dimension];
        else return value;
    }
    inline double get_volume()
    {
        double V=1;
        for (unsigned int i = 0; i < this->boundaries.size(); ++i)
        {
            V*=this->boundaries[i];
        }
        return V;
    }
    std::vector<double> boundaries;
};

#endif
