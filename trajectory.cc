#include "trajectory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


Trajectory::Trajectory()
{

}

//! read a LAMMPS atom file 
void Trajectory::read_atom(string path)
{
    constexpr unsigned int d = 3;
    constexpr unsigned int ATOM_HEADER_SIZE = 6+d;
    
    Configuration* ref_config = nullptr;
    ifstream in(path);
    
    int countlines=0;
    string line;

    while(in.good()){
        getline(in,line);
        countlines++;
    }
    in.close();

    // reopen
    in.open(path);

    // Get any other configurations in this file (i.e. in the trajectory).
    this->sequence.push_back( Configuration() );
    int counter=0;
    while (in)
    {
        if (ref_config) this->sequence[counter].read_atom(in, *ref_config);
        else
        {   
            // Do a blind read, assigning the number of particles, the species etc.
            this->sequence[counter].read_atom(in);
            
            this->num_particles = this->sequence[counter].system_size();
            unsigned int num_frames = countlines/(this->num_particles+ATOM_HEADER_SIZE);
            this->sequence.resize(num_frames);
            ref_config = &this->sequence[counter];
        }
        
        // Get the next character to trigger the eof flag if we're at the end.
        counter++;
        in.get();
    }
    
    in.close();
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

void Trajectory::print_configuration(unsigned int frame)
{
    cerr << frame << endl;
    // this->sequence[frame].print_neighbours();
}


void Trajectory::compute_neighbour_correlation(bool sorting)
{
    cerr << sorting << endl;
    // TO DO: TO SPEED UP USING ITERATORS!!!!!!!
   //  for (unsigned int t = 0; t < this->sequence.size()-1; ++t)
   //  {   
   //      for (unsigned int tt = t+1; tt < this->sequence.size(); ++tt)
   //      {
            
   //          this->neigh_corr[tt-t]+=sequence[t].neighbour_overlap(sequence[tt],sorting);
   //          this->neigh_norm[tt-t]++;
   //      }
   //  }
   //  for (unsigned int i = 0; i < neigh_corr.size(); ++i)
   //  {
   //      this->neigh_corr[i]/=neigh_norm[i];
   //  }
   // this->neigh_corr[0]=1;
 }

void Trajectory::save_neighbour_correlation(std::string filename){
    std::ofstream fout(filename);

    for (unsigned int i = 0; i < this->sequence.size(); ++i)
    {
        fout << i << '\t' << this->neigh_corr[i] << "\t" << this->neigh_norm[i] << "\n";
    }
    fout.close();
}

void Trajectory::compute_msd_isf(double q)
{
    constexpr unsigned int d = 3;
    
    this->isf.resize(this->sequence_length());
    this->num_samples.resize(this->sequence_length());
    this->msd.resize(this->sequence_length());
    
    for (unsigned int t = 0; t < this->sequence_length(); ++t) (msd[t]).resize(d);
    
    for (unsigned int t = 0; t < this->sequence_length(); ++t)
    {
        this->isf[t]=0;
        this->num_samples[t]=0;
        for (unsigned int c = 0; c < d; ++c) this->msd[t][c]=0;
    }
    
    vector< vector<double> > msd_isf_table( this->sequence_length(), vector<double>(d+1) );
    for (unsigned int t = 0; t < this->sequence_length()-1; ++t)
    {
        for (unsigned int tt = t+1; tt < this->sequence_length(); ++tt)
        {
            this->sequence[t].cumulative_msd_isf(msd_isf_table[tt-t], sequence[tt], q);
            this->num_samples[tt-t]++;
        }
    }
    
    // normalise
    for (unsigned int t=0; t < this->sequence_length(); ++t)
    {
        if (num_samples[t]>0)
        {
            for (unsigned int c = 0; c < d; ++c) this->msd[t][c] = msd_isf_table[t][c]/num_samples[t];
            this->isf[t] = msd_isf_table[t][d]/num_samples[t];
      }
    }
}

void Trajectory::save_msd_isf(string path)
{
    constexpr int d = 3;
    
    ofstream fout(path);
    // skip the 0th lag-time
    for (unsigned int i = 1; i < this->sequence_length(); ++i)
    {
        fout << i << '\t';
        for (unsigned int c = 0; c < d; ++c) fout << this->msd[i][c] << "\t";
        fout << this->isf[i] << "\t" << num_samples[i] << "\n";
    }
    
    fout.close();
}

void Trajectory::compute_g(unsigned int num_bins, double delta_r)
{
    this->g.resize(num_bins);
    this->delta_bin=delta_r;
    for (auto t=this->sequence.begin(); t!=this->sequence.end(); ++t)
        (*t).cumulative_radial_distribution(this->g, delta_r);
}

void Trajectory::save_g(std::string filename)
{
    std::ofstream fout(filename);
    
    for (unsigned int bin = 0; bin < this->g.size(); ++bin)
    {
        g[bin] /= this->sequence_length();
        fout << (bin+0.5)*this->delta_bin << "\t" << g[bin] << "\n";
    }
    
    fout.close();
}

