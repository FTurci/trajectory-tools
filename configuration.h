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
    void read_xyz(std::istream& in);
    //void read_xyzr(std::string path);
    //void read_pdb(std::string path);
    //void read_lammps(std::string path);
    // Optimised reading utilities for when the number of particles/species are known in advance, thus avoiding potentially expensive dynamic memory allocation.
    void read_xyz(std::string path, const std::vector<unsigned int>& species_distribution);
    void read_xyz(std::istream& in, const std::vector<unsigned int>& species_distribution);
    
    //
    inline const std::vector<unsigned int>& get_dispersity() const
    {
        return this->dispersity;
    }
    
    void read_neighbours(std::string path);
    // apply periodic boundaries to values
    void periodic_boundaries(std::vector<double> &values);
    
    void print_positions(std::ostream& out) const;
    void print_neighbours(int first_particle, int last_particle) const;
    // print all
    void print_neighbours() const;
    
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
    const std::vector<double>& radial_distribution(unsigned int num_bins, double bin_width);
    
protected:
    // Particle data.
    unsigned int num_particles;
    std::vector< Species<3> > particles;
    // A summary of the size of the population of each species in the above vector.
    // NB: These values should essentially be particles[i].size(), but we keep them separate for convenience and to quickly return the size of each population without recalculation.
    std::vector<unsigned int> dispersity;
    
    // Bookkeeping: each particle is assigned a unique id so we can keep track of them individually if need be.
    struct ParticleIndex
    {
        unsigned int species, index;
    };
    std::vector<ParticleIndex> particle_table;
    
    // ...
    std::vector< std::vector<int> > neighbour_table;
    // g(r)
    std::vector<double> g;
    double g_bin_width;
    // experimental g(r)
    std::vector<double> experimental_g;
    
    // number of particles
    //int Npart;
    // number density
    //double density;
    // box sizes
    std::vector<double> box;
};

#endif

