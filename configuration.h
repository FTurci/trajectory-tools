#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <iostream>
#include <vector>
#include <string>

#include "species.h"


// Only supports canonical configurations for now, i.e. particle number N is fixed.
class Configuration
{
public:
    Configuration();
    
    // 'Blind' read functions, which read position data without foreknowledge of number of particles/species etc.
    // NB: these will be a lot slower than preallocating numbers of particles in other read functions, so should only be used on the first configuration within a trajectory.
    void read_xyz(std::string path);
    //void read_xyzr(std::string path);
    //void read_pdb(std::string path);
    //void read_lammps(std::string path);
    // Optimised reading utilities for when the number of particles/species are known in advance, thus avoiding potentially expensive dynamic memory allocation.
    void read_xyz(std::string path, const std::vector<unsigned int>& species_distribution);
    
    //
    inline const std::vector<unsigned int>& get_dispersity() const
    {
        return this->dispersity;
    }
    
    void read_neighbours(std::string file);
    // apply periodic boundaries to values
    void periodic_boundaries(std::vector<double> &values);
    
    void print_positions(std::ostream& out) const;
    void print_neighbours(int first_particle, int last_particle);
    // print all
    void print_neighbours();
    
    // Default stream overloads work on the positions.
    inline friend std::ostream& operator<< (std::ostream& out, const Configuration& config)
    {
        config.print_positions(out);
        return out;
    }
    
    // Pair correlations:
    // compute the average overlap between the lists of neighbours
    double neighbour_overlap(Configuration b, bool sorting=false);
    // exmperimental g(r)
    void experimental_radial_distr();
    // simulation g(r)
    void radial_distr(int nbins, double biwidth);
    
protected:
    unsigned int numParticles;
    std::vector< Species<3> > particles;
    std::vector<unsigned int> dispersity;
    
    std::vector< std::vector<int> > neighbour_table;
    // g(r)
    std::vector<double> g;
    // experimental g(r)
    std::vector<double> experimental_g;
    
    // number of particles
    //int Npart;
    // number density
    double density;
    // box sizes
    std::vector<double> box;
    
    std::vector<double> coordinates;
};

#endif

