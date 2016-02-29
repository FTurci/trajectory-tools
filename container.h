#ifndef __CONTAINER_H
#define __CONTAINER_H

#include <array>


class Container
{
public:
    const static int d = 3;

    Container() = default;
    Container(const Container&) = default;
    Container(Container&&) = default;
    Container& operator=(const Container&) = default;
    Container& operator=(Container&&) = default;

    inline std::array<double,d> get_size() const
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
    inline std::array<double,d> get_origin() const
    {
        return this->origin;
    }
    inline void set_origin(std::array<double,d> new_origin)
    {
        for (int c = 0; c < d; ++c) this->origin[c] = new_origin[c] - this->origin[c];
    }
    inline void set_origin(double* new_origin)
    {
        for (int c = 0; c < d; ++c) this->origin[c] = new_origin[c] - this->origin[c];
    }

    inline double apply_boundaries(double value, unsigned int dimension) const
    {
        if (value > this->boundaries[dimension]*0.5) return value-this->boundaries[dimension];
        else if (value < -this->boundaries[dimension]*0.5) return value+this->boundaries[dimension];
        else return value;
    }

protected:
    std::array<double,d> origin;
    std::array<double,d> boundaries;
};

#endif
