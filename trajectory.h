#ifndef __TRAJECTORY_H
#define __TRAJECTORY_H

#include "configuration.h"


class trajectory
{
public:
    trajectory();
    // void read_sequence(std::vector<string> config_paths, std::vector<string> neighbour_paths=std::vector<string>());
    void read_sequence_neighbours(std::vector<std::string> neighbour_paths);
    void print_configuration(int frame);

    void compute_neighbour_correlation(bool sorting);
    void save_neighbour_correlation(std::string);
    int length();
private:

    std::vector<configuration> sequence;
    std::vector<double> neigh_corr;
    std::vector<double> neigh_norm;

    /* data */
};

#endif

