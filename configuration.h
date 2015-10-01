#ifndef __configuration_H
#define __configuration_H


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>  
#include <algorithm>



class configuration
{
public:
    configuration();

    void read_neighbours(std::string file);
    void print_neighbours(int first_particle,int last_particle);
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

