#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <string>

#include "species.h"
#include "container.h"


//class RadialDistribution

// Only supports canonical configurations for now, i.e. particle number N is fixed.
class Configuration : public Container
{
public:
    constexpr static int d = 3;

    Configuration();
    Configuration(const Configuration&) = default;
    Configuration(Configuration&&) = default;
    Configuration& operator=(const Configuration&) = default;
    Configuration& operator=(Configuration&&) = default;

    // Summary information (number of particles, species breakdown etc).
    struct Diversity
    {
        Diversity() : system_size(0) { }
        Diversity(const Diversity&) = default;
        Diversity(Diversity&&) = default;
        Diversity& operator=(const Diversity&) = default;
        Diversity& operator=(Diversity&&) = default;

        unsigned int system_size;
        // The name and population of each species.
        std::vector< std::pair<std::string,unsigned int> > species;
        // A map from the above name -> internal index
        std::map<std::string, unsigned int> species_map;

        inline friend std::ostream& operator<< (std::ostream& out, const Diversity& diversity)
        {
            out << "total: " << diversity.system_size << "\n";
            for (unsigned int i = 0; i < diversity.species.size(); ++i)
                out << diversity.species[i].first << ": " << diversity.species[i].second << "\n";
            return out;
        }
    };
    inline const Diversity& summary() const
    {
        return this->diversity;
    }

    // Particle access.
    struct ParticleIndex
    {
        unsigned int species, index;
    };
    inline const ParticleIndex& operator[] (int n) const
    {
        return this->particle_table[n];
    }
    inline const ParticleIndex& operator() (int n) const
    {
        return this->particle_table[n];
    }
    inline const double* operator() (int s, int n) const
    {
        return &this->particles[s][n];
    }

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
    void read_xyz(std::string path, const Configuration& ref_config);
    void read_xyz(std::istream& in, const Configuration& ref_config);
    void read_xyz(std::string path, const Diversity& ref_div);
    void read_xyz(std::istream& in, const Diversity& ref_div);
    void read_atom(std::string path, const Configuration& ref_config);
    void read_atom(std::istream& in, const Configuration& ref_config);
    //void read_atom(std::string path, const std::vector<unsigned int>& species_distribution);
    //void read_atom(std::istream& in, const std::vector<unsigned int>& species_distribution);

    // Write functions.
    void write_xyz(std::string path) const;
    void write_xyz(std::ostream& out) const;

    // ************************ OUT-OF-DATE
    inline unsigned int system_size() const
    {
        return this->num_particles;
    }
    inline const std::vector<unsigned int>& get_dispersity() const
    {
        return this->dispersity;
    }
    // ************************ TO BE DELETED

    Configuration subsystem(double distance) const;
    Configuration subsystem(double distance, unsigned int central_particle) const;
    Configuration subsystem(double distance, const double* center) const;
    Configuration exclude_subsystem(double distance) const;
    Configuration exclude_subsystem(double distance, unsigned int central_particle) const;
    Configuration exclude_subsystem(double distance, const double* center) const;
    void set_origin(std::array<double,d> new_origin);

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

    // Computes minimal box needed to enclose all the particles, useful when the exact box is not known (when e.g. reading xyz files).
    // Note: accuracy improves with larger/denser systems.
    std::vector<double> pseudo_box() const;

    // Pair correlations:
    // compute the average overlap between the lists of neighbours
    double neighbour_overlap(Configuration& b, bool sorting=false);

    std::vector<double> msd_isf(const Configuration& b, const double q) const;
    void cumulative_msd_isf(std::vector<double>& msd_isf_total, const Configuration& b, const double q) const;
    // exmperimental g(r)
    void experimental_radial_distr() const;
    // Simulation g(r).
    std::vector<double> radial_distribution(unsigned int num_bins, double bin_width) const;
    std::vector<double> radial_distribution(unsigned int species, unsigned int num_bins, double bin_width) const;
    std::vector<double> radial_distribution(unsigned int species_a, unsigned int species_b, unsigned int num_bins, double bin_width) const;
    // Cumulative calculations for computing averages.
    void cumulative_radial_distribution(std::vector<double>& g_total, double bin_width) const;
    void cumulative_radial_distribution(unsigned int species, std::vector<double>& g_total, double bin_width) const;
    void cumulative_radial_distribution(unsigned int species_a, unsigned int species_b, std::vector<double>& g_total, double bin_width) const;

protected:
    // Particle data.
    unsigned int num_particles;
    std::vector< Species<d> > particles;
    // A summary of the size of the population of each species in the above vector.
    // NB: These values should essentially be particles[i].size(), but we keep them separate for convenience and to quickly return the size of each population without recalculation.
    Diversity diversity;
    std::vector<unsigned int> dispersity;

    // Bookkeeping: each particle is assigned a unique id so we can keep track of them individually if need be.
    std::vector<ParticleIndex> particle_table;

    // ...
    std::vector< std::vector<int> > neighbour_table;

    // number density
    //double density;
};

#endif
