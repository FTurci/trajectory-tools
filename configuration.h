#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <vector>
#include <string>

#include "species.h"


// Only supports canonical configurations for now, i.e. particle number N is fixed.
class Configuration
{
public:
    Configuration();
    
    // 'Blind' read function without foreknowledge of number of particles/species etc.
    // NB: this will be a lot slower than preallocating numbers of particles in other read functions, so should only be used on the first configuration within a trajectory.
    void read_xyz(std::string path);
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
    
protected:
    unsigned int numParticles;
    std::vector< Species<3> > particles;
    
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

