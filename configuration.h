#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <vector>
#include <string>


template <unsigned int d>
class Species
{
public:
    Species(unsigned long N) : coords(d*N), N(N) { }
    Species(const Species& copy) : coords(copy.coords), N(copy.N) { }
    
    inline double& operator[] (int n)
    {
        return coords[n];
    }
    inline double operator[] (int n) const
    {
        return coords[n];
    }
    inline double& operator() (int n, int c)
    {
        return coords[n*N+c];
    }
    inline double operator() (int n, int c) const
    {
        return coords[n*N+c];
    }
    
protected:
    std::vector<double> coords;
    const unsigned long N;
};

class Configuration
{
public:
    Configuration();
    
    void read_neighbours(std::string file);
    void print_neighbours(int first_particle, int last_particle);
    // print all
    void print_neighbours();
    // apply periodic boundaries to values
    void periodic_boundaries(std::vector<double> &values);
    
    // Pair correlations:
    // compute the average overlap between the lists of neighbours
    double neighbour_overlap(Configuration b, bool sorting=false);
    // exmperimental g(r)
    void experimental_radial_distr();
    // simulation g(r)
    void radial_distr(int nbins,double biwidth);

private:
    std::vector< std::vector<int> > neighbour_table;
    // g(r)
    std::vector<double> g;
    // experimental g(r)
    std::vector<double> experimental_g;


    // number of particles
    int Npart;
    // number density
    double density;
    // box sizes
    std::vector<double> box;

    std::vector<double> coordinates;

};

#endif

