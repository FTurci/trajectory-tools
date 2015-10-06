#ifndef __TRAJECTORY_H
#define __TRAJECTORY_H

#include "configuration.h"
#include <list>
#include <string>
#include <cmath>

class Trajectory
{
public:
    Trajectory();
   
    void read_atom(std::string path);
    void read_sequence(std::vector<std::string> config_paths,
                       std::vector<std::string> neighbour_paths=std::vector<std::string>());
    void read_sequence_neighbours(std::vector<std::string> neighbour_paths);
    void print_configuration(int frame);

    void compute_msd_isf(double q);
    void save_msd_isf(std::string);

    void compute_g(int num_bins, double delta_r);
    void save_g(std::string);
    
    void compute_neighbour_correlation(bool sorting);
    void save_neighbour_correlation(std::string);
    

    inline int length(){return this->sequence.size();};


private:
    std::vector<Configuration> sequence;
    std::vector<double> neigh_corr;
    std::vector<double> neigh_norm;
    std::vector<double> isf;
    std::vector<std::vector<double> > msd;
    std::vector<int> num_samples;
    std::vector<double> g;

    int num_particles;
    double delta_bin;

    /* data */
};

#endif

