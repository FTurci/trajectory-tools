#ifndef __CONFIGURATION_BASE_H
#define __CONFIGURATION_BASE_H

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <string>

#include "species.h"
#include "container.h"
#include "subsystems.h"


// Only supports canonical configurations for now, i.e. particle number N is fixed.
class BaseConfiguration : public Container
{
public:
    BaseConfiguration() = default;
    BaseConfiguration(const Container&);
    BaseConfiguration(const BaseConfiguration&) = default;
    BaseConfiguration(BaseConfiguration&&) = default;
    BaseConfiguration& operator=(const BaseConfiguration&) = default;
    BaseConfiguration& operator=(BaseConfiguration&&) = default;

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
    inline unsigned int size() const
    {
        return this->diversity.system_size;
    }
    inline unsigned int num_species() const
    {
        return this->diversity.species.size();
    }
    inline bool is_empty() const
    {
        return (this->diversity.system_size == 0);
    }

    virtual void read_xyz(std::string path);
    virtual void read_xyz(std::istream& in);

    // Write functions.
    void write_coords(std::string path) const;
    void write_coords(std::ostream& out) const;
    void write_xyz(std::string path) const;
    void write_xyz(std::ostream& out) const;

    // Default stream overloads assume xyz files.
    inline friend std::ostream& operator<< (std::ostream& out, const BaseConfiguration& config)
    {
        config.write_xyz(out);
        return out;
    }

    //void periodic_boundaries(std::vector<double> &values);

    // Computes minimal box needed to enclose all the particles, useful when the exact box is not known (when e.g. reading xyz files).
    // Note: accuracy improves with larger/denser systems.
    std::array<double,d> pseudo_box() const;

    // Pair correlations for the 3-sphere
    /* Pair correlations:
    void read_neighbours(std::string path);
    // compute the average overlap between the lists of neighbours
    double neighbour_overlap(Configuration& b, bool sorting=false);
    // printing
    void print_neighbours(int first_particle, int last_particle) const;
    void print_neighbours() const;*/

    /*std::vector<double> msd_isf(const Configuration& b, const double q) const;
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
    void cumulative_radial_distribution(unsigned int species_a, unsigned int species_b, std::vector<double>& g_total, double bin_width) const;*/

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
        return position(s,n);
    }
    
    virtual const double* position(int s, int n) const = 0;

    // Generic operations on particles using functors.

    template<class Operation>
    void for_each(Operation& operation) const;
    template<class Operation>
    void for_each_pair(Operation& operation) const;

protected:
    virtual double* position(int s, int n) = 0;
    virtual inline double* operator() (int s, int n)
    {
        return position(s,n);
    }

    // Bookkeeping: each particle is assigned a unique id, diversity tracks
    // the number of particles the object can access in the Species vector
    // (NB: this may be a sub-configuration).
    Diversity diversity;
    std::vector<ParticleIndex> particle_table;

    // ...
    //std::vector< std::vector<int> > neighbour_table;

    // number density
    //double density;
    template<class View, class Criteria> View subsystem(Criteria criteria) const;
    //template<class T> T exclude_subsystem(double distance) const
    //{
    //}
};

template<class Operation>
void BaseConfiguration::for_each(Operation& operation) const
{
    const ParticleIndex* lookup;
    const double* r;

    for (unsigned int i = 0; i < this->size()-1; ++i)
    {
        lookup = &this->particle_table[i];
        r = this->position(lookup->species, lookup->index);
        operation(r);
    }
}

template<class Operation>
void BaseConfiguration::for_each_pair(Operation& operation) const
{
    const ParticleIndex* lookup;
    const double *r1, *r2;

    for (unsigned int i = 0; i < this->size()-1; ++i)
    {
        lookup = &this->particle_table[i];
        r1 = this->position(lookup->species, lookup->index);

        for (unsigned int j = i+1; j < this->size(); ++j)
        {
            lookup = &this->particle_table[j];
            r2 = this->position(lookup->species, lookup->index);
            operation(r1, r2);
        }
    }
}

template<class View, class Criteria>
View BaseConfiguration::subsystem(Criteria criteria) const
{
    View sub(this);
    sub.boundaries = this->boundaries;
    sub.diversity.system_size = 0;

    for (unsigned int i = 0; i < this->size(); ++i)
    {
        auto lookup = &this->particle_table[i];
        const double* r = this->position(lookup->species, lookup->index);

        if (criteria(r))
        {
            sub.positions[lookup->species].push_back(const_cast<double*>(r));

            // Update the bookkeeping.
            sub.particle_table.push_back(ParticleIndex());
            ParticleIndex& id = sub.particle_table.back();
            id.species = lookup->species;
            id.index = sub.diversity.species[lookup->species].second;

            sub.diversity.species[lookup->species].second++;
            sub.diversity.system_size++;
        }
    }

    // NB: Will use more memory than necessary while this is commented out,
    // but then there won't be performance overheads from reallocating memory.
    //sub.particle_table.shrink_to_fit();
    //for (unsigned int i = 0; i < this->num_species(); ++i)
    //    this->positions[i].shrink_to_fit();

    return sub;
}


#endif
