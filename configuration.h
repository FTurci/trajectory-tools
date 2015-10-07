#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <iostream>
#include <vector>
#include <string>

#include "species.h"
#include "container.h"


// Only supports canonical configurations for now, i.e. particle number N is fixed.
class Configuration : public Container
{
public:
    Configuration();
    Configuration(const Configuration& copy);
    //Configuration(Configuration&& move) = default;
    //Configuration& operator=(const Configuration& copy) & = default;
    //Configuration& operator=(Configuration&& move) & = default;
    
    // 'Blind' read functions, which read position data without foreknowledge of number of particles/species etc.
    // NB: these will be a lot slower than preallocating numbers of particles in other read functions, so should only be used on the first configuration within a trajectory.
    void read_xyz(std::string path);
    void read_xyz(std::istream& in);
    void read_atom(std::string path);
    void read_atom(std::istream& in);
    //void read_xyzr(std::string path);
    //void read_pdb(std::string path);
    //void read_lammps(std::string path);
    // Optimised reading utilities for when the configuration information is known in advance, to avoid potentially expensive dynamic memory allocation.
    void read_xyz(std::string path, const std::vector<unsigned int>& species_distribution);
    void read_xyz(std::istream& in, const std::vector<unsigned int>& species_distribution);
    void read_atom(std::string path, const Configuration& ref_config);
    void read_atom(std::istream& in, const Configuration& ref_config);
    //void read_atom(std::string path, const std::vector<unsigned int>& species_distribution);
    //void read_atom(std::istream& in, const std::vector<unsigned int>& species_distribution);
    
    //
    inline const std::vector<unsigned int>& get_dispersity() const
    {
        return this->dispersity;
    }
    inline unsigned int system_size()
    {
        return this->num_particles;
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
    double neighbour_overlap(Configuration& b, bool sorting=false);
    
    std::vector<double> msd_isf(const Configuration& b, const double q) const;
    void cumulative_msd_isf(std::vector<double>& msd_isf_total, const Configuration& b, const double q) const;
    // exmperimental g(r)
    void experimental_radial_distr() const;
    // simulation g(r)
    std::vector<double> radial_distribution(unsigned int num_bins, double bin_width) const;
    // Cumulative g(r) for computing averages.
    void cumulative_radial_distribution(std::vector<double>& g_total, double bin_width) const;
    
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
    
    // number density
    //double density;
};

#endif

