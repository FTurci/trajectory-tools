#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include <vector>
#include <string>


class configuration
{
public:
    configuration();

    void read_neighbours(std::string file);
    void print_neighbours(int first_particle, int last_particle);
    // print all
    void print_neighbours();

    double neighbour_overlap(configuration b, bool sorting=false);
    /* data */
private:
    std::vector< std::vector<int> > Neighbour_Table;
    // number of particles
    int Npart;
};

#endif

