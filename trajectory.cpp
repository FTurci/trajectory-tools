#include "trajectory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


trajectory::trajectory()
{

}

void trajectory::read_sequence(std::vector<string> config_paths, std::vector<string> neighbour_paths)
{
}

void trajectory::read_sequence_neighbours(vector<string> path_list)
{
    for (int i = 0; i < path_list.size(); ++i)
    {
        configuration C;
        C.read_neighbours(path_list[i]);
        this->sequence.push_back(C);
    }
    this->neigh_corr.resize(sequence.size());
    this->neigh_norm.resize(sequence.size());
}

void trajectory::print_configuration(int frame)
{
    this->sequence[frame].print_neighbours();
}


void trajectory::compute_neighbour_correlation(bool sorting)
{
    for (int t = 0; t < this->sequence.size()-1; ++t)
    {   
        for (int tt = t+1; tt < this->sequence.size(); ++tt)
        {
            // std::cout<<"T "<<t<<"  TT "<<tt<<std::endl;
            this->neigh_corr[tt-t]+=sequence[t].neighbour_overlap(sequence[tt],sorting);
            this->neigh_norm[tt-t]++;
        }
    }
    for (int i = 0; i < neigh_corr.size(); ++i)
    {
        this->neigh_corr[i]/=neigh_norm[i];
    }
   this->neigh_corr[0]=1;
 }

void trajectory::save_neighbour_correlation(std::string filename){
    std::ofstream fout(filename);

    for (int i = 0; i < this->sequence.size(); ++i)
    {
        fout<<i<<'\t'<<this->neigh_corr[i]<<"\t"<<this->neigh_norm[i]<<std::endl;
    }
    fout.close();
}

int trajectory::length(){
    return this->sequence.size();
}

