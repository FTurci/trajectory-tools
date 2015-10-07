#ifndef __SPECIES_H
#define __SPECIES_H

#include <vector>


template <unsigned int d>
class Species
{
public:
    Species(unsigned long N = 0) : coords(d*N), N(N) { }
    Species(const Species& copy) : coords(copy.coords), N(copy.N) { }
    Species(const std::vector<double>& copy) : coords(copy), N(copy.size()/d) { }

    inline double& operator[] (int n)
    {
        return this->coords[n*d];
    }
    inline const double& operator[] (int n) const
    {
        return this->coords[n*d];
    }
    inline double& operator() (int n, int c)
    {
        return this->coords[n*d+c];
    }
    inline double operator() (int n, int c) const
    {
        return this->coords[n*d+c];
    }
    inline unsigned long size() const
    {
        return this->N;
    }

protected:
    std::vector<double> coords;
    const unsigned long N;
};

#endif
