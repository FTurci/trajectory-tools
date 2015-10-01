#ifndef __trajectory_H
#define __trajectory_H


#include "configuration.h"

class trajectory
{
public:
    trajectory();
    void read_sequence(std::string stem, int first, int last);
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