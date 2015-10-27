#ifndef __TRAJECTORY_H
#define __TRAJECTORY_H

#include "configuration.h"
#include <list>
#include <string>
#include <cmath>

class Trajectory
{
public:
    //constexpr unsigned int d = 3;
    Trajectory();

    //void read_xyz(std::vector<std::string> config_paths);
    void read_atom(std::string path);
    void read_sequence_neighbours(std::vector<std::string> neighbour_paths);
    void print_configuration(unsigned int frame);

    void compute_msd_isf(double q);
    void print_msd_isf(std::string path);
    void print_msd_isf(std::ostream& out);

    void compute_g(unsigned int num_bins, double delta_r);
    void print_g(std::string path);
    void print_g(std::ostream& out);

    void compute_neighbour_correlation(bool sorting);
    void save_neighbour_correlation(std::string);

    inline std::vector<double> container_size() const
    {
        return this->sequence[0].get_size();
    }
    inline unsigned int system_size()
    {
        return this->num_particles;
    }
    inline unsigned int sequence_length()
    {
        return this->sequence.size();
    }

private:
    std::vector<Configuration> sequence;
    std::vector<double> neigh_corr;
    std::vector<double> neigh_norm;
    std::vector<unsigned int> num_samples;

    std::vector<double> isf;
    std::vector<std::vector<double> > msd;
    std::vector<std::vector<double> > g;

    unsigned int num_particles;
    double delta_bin;

    /* data */
};

#endif
