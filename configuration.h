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

    // Pair correlations:
    // compute the average overlap between the lists of neighbours
    double neighbour_overlap(configuration b, bool sorting=false);
    // exmperimental g(r)
    void experimental_pair_correlation();
    // simulation g(r)
    void pair_correlation();

private:
    std::vector< std::vector<int> > neighbour_table;
    // g(r)
    std::vector<double> g;
    // experimental g(r)
    std::vector<double> experimental_g;


    // number of particles
    int Npart;
    // box sizes
    std::vector<double> box;

    std::vector<double> coordinates;

};

#endif

