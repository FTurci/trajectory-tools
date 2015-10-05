#include "trajectory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


Trajectory::Trajectory()
{

}

void Trajectory::read_sequence(std::vector<string> config_paths, std::vector<string> neighbour_paths)
{
    for (auto it = config_paths.begin(); it != config_paths.end(); ++it) cout << *it << "\n";
    for (auto it = neighbour_paths.begin(); it != neighbour_paths.end(); ++it) cout << *it << "\n";
}

void Trajectory::read_sequence_neighbours(vector<string> path_list)
{
    // Error checking: this sequence must be compatible with any sequences already loaded, otherwise it cannot be part of the same trajectory.
    if (this->sequence.size())
    {
    }
    
    for (unsigned int i = 0; i < path_list.size(); ++i)
    {
        Configuration C;
        C.read_neighbours(path_list[i]);
        this->sequence.push_back(C);
    }
    this->neigh_corr.resize(sequence.size());
    this->neigh_norm.resize(sequence.size());
}

void Trajectory::print_configuration(int frame)
{
    this->sequence[frame].print_neighbours();
}


void Trajectory::compute_neighbour_correlation(bool sorting)
{
    for (unsigned int t = 0; t < this->sequence.size()-1; ++t)
    {   
        for (unsigned int tt = t+1; tt < this->sequence.size(); ++tt)
        {
            // cout<<"T "<<t<<"  TT "<<tt<< "\n";
            this->neigh_corr[tt-t]+=sequence[t].neighbour_overlap(sequence[tt],sorting);
            this->neigh_norm[tt-t]++;
        }
    }
    for (unsigned int i = 0; i < neigh_corr.size(); ++i)
    {
        this->neigh_corr[i]/=neigh_norm[i];
    }
   this->neigh_corr[0]=1;
 }

void Trajectory::save_neighbour_correlation(std::string filename){
    std::ofstream fout(filename);

    for (unsigned int i = 0; i < this->sequence.size(); ++i)
    {
        fout<<i<<'\t'<<this->neigh_corr[i]<<"\t"<<this->neigh_norm[i]<<"\n";
    }
    fout.close();
}

int Trajectory::length(){
    return this->sequence.size();
}

