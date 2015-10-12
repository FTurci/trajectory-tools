#include "trajectory.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "utilities.h"


Trajectory::Trajectory()
{

}

//! read a LAMMPS atom file
void Trajectory::read_atom(string path)
{
    constexpr unsigned int d = 3;
    constexpr unsigned int ATOM_HEADER_SIZE = 6+d;

    // Count how many lines there are, so after reading the first configuration we know how many frames there are in this trajectory.
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open file ", path);
    unsigned long num_lines = 0;
    string line;
    while(in)
    {
        getline(in, line);
        num_lines++;
    }
    in.close();

    // Open for reading this time.
    in.open(path);

    // Get any other configurations in this file (i.e. in the trajectory).
    Configuration ref_config;
    unsigned int counter = 0;
    while (in)
    {
        if (counter) this->sequence[counter].read_atom(in, ref_config);
        else
        {
            // Do a blind read, assigning the number of particles, the species etc.
            ref_config.read_atom(in);
            this->num_particles = ref_config.system_size();
            // Create space for the rest of the trajectory.
            unsigned int num_frames = num_lines/(this->num_particles+ATOM_HEADER_SIZE);
            this->sequence = vector<Configuration>(num_frames);
            this->sequence[counter] = ref_config;
        }

        // Get the next character to trigger the eof flag if we're at the end.
        counter++;
        in.get();
    }

    in.close();
}

/*void Trajectory::read_sequence(std::vector<string> config_paths, std::vector<string> neighbour_paths)
{
    for (auto it = config_paths.begin(); it != config_paths.end(); ++it) cout << *it << "\n";
    for (auto it = neighbour_paths.begin(); it != neighbour_paths.end(); ++it) cout << *it << "\n";
    }*/

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
            this->num_samples[tt-t]+=this->num_particles;
        }
    }

    // normalise
    for (unsigned int t = 0; t < this->sequence_length(); ++t)
    {
        if (num_samples[t] > 0)
        {
            for (unsigned int c = 0; c < d; ++c) this->msd[t][c] = msd_isf_table[t][c]/num_samples[t];
            this->isf[t] = msd_isf_table[t][d]/num_samples[t];
      }
    }
}

void Trajectory::print_msd_isf(string path)
{
    ofstream out(path);
    if (!out) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path, " for writing");
    this->print_msd_isf(out);
    out.close();
}

void Trajectory::print_msd_isf(ostream& out)
{
    constexpr int d = 3;

    // skip the 0th lag-time
    for (unsigned int i = 1; i < this->sequence_length(); ++i)
    {
        out << i << '\t';
        for (unsigned int c = 0; c < d; ++c) out << this->msd[i][c] << "\t";
        out << this->isf[i] << "\t" << num_samples[i] << "\n";
    }
}

void Trajectory::compute_g(unsigned int num_bins, double delta_r)
{
    const unsigned int num_species = this->sequence[0].get_dispersity().size();
    // We need to make room for the g(r) for each species and the total.
    this->g = vector< vector<double> >(num_species+1, vector<double>(num_bins) );
    this->delta_bin = delta_r;

    for (auto t = this->sequence.begin(); t != this->sequence.end(); ++t)
    {
        // g(r) for individual species.
        for (unsigned int species = 0; species < num_species; ++species)
            (*t).cumulative_radial_distribution(species, this->g[species], delta_r);
        // Add the missing terms to the combined form: We can add the individual contributions calculated above later.
        for (unsigned int species_a = 0; species_a < num_species-1; ++species_a)
            for (unsigned int species_b = species_a+1; species_b < num_species; ++species_b)
                (*t).cumulative_radial_distribution(species_a, species_b, this->g[num_species], delta_r);
    }

    // Fill in the missing information for the combined form.
    for (unsigned int species = 0; species < num_species; ++species)
        for (unsigned int bin = 0; bin < num_bins; ++bin)
            this->g[num_species][bin] += this->g[species][bin];

    // We want the average so we have to divide by the number of samples.
    for (unsigned int species = 0; species <= num_species; ++species)
        for (unsigned int bin = 0; bin < num_bins; ++bin)
            this->g[species][bin] /= this->sequence_length();
}

void Trajectory::print_g(string path)
{
    ofstream out(path);
    if (!out) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path, " for writing");
    this->print_g(out);
    out.close();
}

void Trajectory::print_g(ostream& out)
{
    for (unsigned int bin = 0; bin < this->g[0].size(); ++bin)
    {
        out << (bin+0.5)*this->delta_bin;
        for (auto species = this->g.begin(); species != this->g.end(); ++species)
            out << "\t" << (*species)[bin];
        out << "\n";
    }
}
