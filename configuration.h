#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include "configuration_base.h"
#include "configuration_view.h"


// Only supports canonical configurations for now, i.e. particle number N is fixed.
class Configuration : public BaseConfiguration
{
public:
    Configuration();
    Configuration(const Configuration&) = default;
    Configuration(Configuration&&) = default;
    Configuration& operator=(const Configuration&) = default;
    Configuration& operator=(Configuration&&) = default;

    // 'Blind' read functions, which read position data without foreknowledge of number of particles/species etc.
    // NB: these will be slower than preallocating numbers of particles in other read functions, so should only be used on the first configuration within a trajectory.
    void read_xyz(std::string path);
    void read_xyz(std::istream& in);
    // Optimised reading utilities for when the configuration information is known in advance, to avoid potentially expensive dynamic memory allocation.
    void read_xyz(std::string path, const Diversity& ref_div);
    void read_xyz(std::string path, const BaseConfiguration* ref_conf);
    void read_xyz(std::istream& in, const Diversity& ref_div);
    void read_xyz(std::istream& in, const BaseConfiguration* ref_conf);

    // Same for LAMMPS atom files.
    void read_atom(std::string path);
    void read_atom(std::istream& in);
  /*void read_atom(std::string path, const Configuration& ref_config);
    void read_atom(std::istream& in, const Configuration& ref_config);
    void read_atom(std::string path, const std::vector<unsigned int>& species_distribution);
    void read_atom(std::istream& in, const std::vector<unsigned int>& species_distribution);*/

    //void read_xyzr(std::string path);
    //void read_pdb(std::string path);
    //void read_lammps(std::string path);

    //template<class View> View BaseConfiguration::subsystem(double distance, const double* from) const;
    template<class Criteria> ConfigurationView subsystem(Criteria criteria) const
    {
        return BaseConfiguration::subsystem<ConfigurationView>(criteria);
    }
    //Configuration subsystem(double distance, unsigned int central_particle) const;
    //Configuration subsystem(double distance, const double* center) const;
    //Configuration exclude_subsystem(double distance) const;
    //Configuration exclude_subsystem(double distance, unsigned int central_particle) const;
    //Configuration exclude_subsystem(double distance, const double* center) const;
    //void set_origin(std::array<double,d> new_origin);

    inline const double* position(int s, int n) const
    {
        return &(this->particles)[s][n];
    }
    inline double* position(int s, int n)
    {
        return &(this->particles)[s][n];
    }
protected:
    // Particles and their access.
    std::vector< Species<d> > particles;
};

#endif
