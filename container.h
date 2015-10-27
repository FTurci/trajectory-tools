#ifndef __CONTAINER_H
#define __CONTAINER_H

#include <vector>


class Container
{
public:
    Container();
    Container(const Container& copy);

    inline std::vector<double> get_size() const
    {
        return this->boundaries;
    }

    inline double get_volume() const
    {
        double V=1;
        for (unsigned int i = 0; i < this->boundaries.size(); ++i)
        {
            V *= this->boundaries[i];
        }
        return V;
    }

protected:
    inline double apply_boundaries(double value, unsigned int dimension) const
    {
        if(value > this->boundaries[dimension]*0.5) return value-this->boundaries[dimension];
        else if(value < -this->boundaries[dimension]*0.5) return value+this->boundaries[dimension];
        else return value;
    }
    std::vector<double> boundaries;
};

#endif
