#include <iostream>
#include <fstream>
#include <sstream>  
#include <algorithm>
#include <map>
using namespace std;

#include "configuration.h"
#include "utilities.h"


Configuration::Configuration() : num_particles(0), g_bin_width(0.)
{
}

void Configuration::read_xyz(string path)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_xyz(in);
    in.close();
}

void Configuration::read_xyz(istream& in)
{
    const int d = 3;
    
    if (this->particles.size()) throw Exception(__PRETTY_FUNCTION__, ": attempting to load new configuration into a preexisting configuration");
        
    // The first line states the number of particles.
    string line;
    getline(in, line);
    this->num_particles = stoi(line);
    this->particle_table = vector<ParticleIndex>(this->num_particles);
    
    // The comment line can just be discarded (for now, this may change in later revisions).
    getline(in, line);
    
    // Declare all variables outside of the loop for optimisation.
    // NB: this could be unnecessary as this is not meant to be a fast function.
    string species;
    int species_index;
    map<string, int> index_list;
    map<string, int>::iterator found;
    vector< vector<double> > positions;
    double x[d];
    for (unsigned int n = 0; n < this->num_particles; ++n)
    {
        in >> species;
        // If this is a new species then we have to add new data structures for it.
        found = index_list.find(species);
        if (found == index_list.end())
        {
            positions.push_back( vector<double>() );
            species_index = index_list.size();
            index_list[species] = species_index;
        }
        else species_index = found->second;
        
        // Get the coordinates and add them to the list.
        for (unsigned int c = 0; c < d; c++) in >> x[c];
        positions[species_index].insert(positions[species_index].end(), x, x+d);
        
        // Bookkeeping so every particle has a unique id (based on their order in the xyz file).
        this->particle_table[n].species = species_index;
        this->particle_table[n].index = (positions[species_index].size()/d)-1;
    }
    
    // Save the coordinate data as static-length species, using the minimum storage space.
    for (auto it = index_list.begin(); it != index_list.end(); ++it)
    {
        species_index = it->second;
        positions[species_index].resize( positions[species_index].size() );
        this->particles.push_back( Species<d>(positions[species_index]) );
        this->dispersity.push_back( positions[species_index].size()/d );
    }
    
    /*************** DEBUG ****************/
    ParticleIndex* id;
    id = &this->particle_table[0];
    this->boundaries = vector<double>(d);
    for (unsigned int c = 0; c < d; ++c) this->boundaries[c] = abs(this->particles[id->species](id->index,c));
    for (unsigned int n = 1; n < this->num_particles; ++n)
    {
        id = &this->particle_table[n];
        for (unsigned int c = 0; c < d; ++c) this->boundaries[c] = max(this->boundaries[c], abs(this->particles[id->species](id->index,c)));
    }
    for (unsigned int c = 0; c < d; ++c) this->boundaries[c] = this->boundaries[c]*2;
}

void Configuration::read_xyz(string path, const vector<unsigned int>& species_distribution)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_xyz(in, species_distribution);
    in.close();
}

void Configuration::read_xyz(std::istream& in, const vector<unsigned int>& species_distribution)
{
    const int d = 3;
    
    if (this->particles.size())
        throw Exception(__PRETTY_FUNCTION__,
                        ": attempting to load new configuration into a preexisting configuration");
    
    // Make sure the given distribution is compatible with any other data structures we have.
    unsigned int n = 0;
    for (auto it = species_distribution.begin(); it != species_distribution.end(); ++it) n += *it;
    if (!this->num_particles)
    {
        this->num_particles = n;
        this->particle_table = vector<ParticleIndex>(this->num_particles);
    }
    else
    {
        if (this->num_particles != n)
            throw Exception(__PRETTY_FUNCTION__, ": attempting to load incompatible configuration: ",
                            "new size=", n, " but expected size=", this->num_particles);
    }
    if (!this->dispersity.empty())
    {
        if (this->dispersity.size() != species_distribution.size())
            throw Exception(__PRETTY_FUNCTION__, ": attempting to load incompatible configuration: ",
                            "new species=", species_distribution.size(), " but expected species=", this->dispersity.size());
        else
        {
            for (unsigned int i = 0; i < this->dispersity.size(); i++)
                if (this->dispersity[i] != species_distribution[i])
                    throw Exception(__PRETTY_FUNCTION__, ": attempting to load incompatible configuration: ",
                                    "different dispersities");
        }
    }
    else this->dispersity = vector<unsigned int>(species_distribution);
    
    // The first line states the number of particles.
    string line;
    getline(in, line);
    if (stoul(line) != this->num_particles)
        throw Exception(__PRETTY_FUNCTION__, ": attempting to read configuration with ", stoi(line),
                        " particles when configuration expects ", this->num_particles);
    
    // The comment line can just be discarded (for now, this may change in later revisions).
    getline(in, line);
    
    // Preallocate the data types.
    for (auto it = this->dispersity.begin(); it != this->dispersity.end(); ++it)
        this->particles.push_back( Species<d>(*it) );
    
    // Declare all variables outside of the loop for optimisation.
    string species;
    int species_index;
    map<string, int> index_list;
    map<string, int>::iterator found;
    vector<unsigned int> count( this->dispersity.size() );
    for (unsigned int n = 0; n < this->num_particles; ++n)
    {
        in >> species;
        // If this is a new species then we have to add new data structures for it.
        found = index_list.find(species);
        if (found == index_list.end())
        {
            species_index = index_list.size();
            index_list[species] = species_index;
            if (index_list.size() > species_distribution.size())
                throw Exception(__PRETTY_FUNCTION__, ": more species found than expected",
                                ": found (at least)=", index_list.size(), ", expected=", species_distribution.size());
        }
        else species_index = found->second;
        
        // Check this new particle does not violate the prescribed dispersity.
        count[species_index] += 1;
        if (count[species_index] > this->dispersity[species_index])
            throw Exception(__PRETTY_FUNCTION__, ": more particles found of species ", species, " than expected",
                            ": found (at least)=", count[species_index], ", expected=", this->dispersity[species_index]);
        
        // Get the coordinates.
        for (unsigned int c = 0; c < d; c++)
            in >> this->particles[species_index]( count[species_index]-1, c );
    }
}

void Configuration::read_neighbours(std::string path)
{
    std::ifstream file(path);
    if(!file.good()) std::cerr << "ERROR: the file " << path << " does not exist\n Forced exit.\n";
    std::vector< std::vector<int> >table;
    for (std::string line; std::getline(file, line); )
    {   
        std::stringstream stream;
        stream.str(line);
        
        int count=0;
        std::vector<int> neighs;
        for (std::string s; std::getline(stream, s,' '); ){
            // skip the first element
            if(count>0) 
            {
                neighs.push_back(StringToNum<int>(s));
            }
            count++;
        }
        // cout<<neighs.size()<<endl;
        this->neighbour_table.push_back(neighs);
    }
    this->num_particles=neighbour_table.size();
}

void Configuration::print_positions(ostream& out) const
{
    const int d = 3;
    for (unsigned int i = 0; i < this->particles.size(); ++i)
    {
        const Species<d>& sp = this->particles[i];
        out << "species " << i << ": " << sp.size() << " particles" << endl;
        for (unsigned int n = 0; n < sp.size(); ++n)
        {
            for (unsigned int c = 0; c < d; ++c)
                out << "  " << sp(n,c);
            out << "\n";
        }
    }
}

void Configuration::print_neighbours(int first, int last) const
{
    for (int i = first; i < last; ++i)
    {
        for (unsigned int j = 0; j < neighbour_table[i].size(); ++j)
            std::cout << this->neighbour_table[i][j] << " ";
        std::cout << "\n";
    }

}

void Configuration::print_neighbours() const
{
    for (unsigned int i = 0; i < this->num_particles; ++i)
    {
        for (unsigned int j = 0; j < neighbour_table[i].size(); ++j)
            std::cout << this->neighbour_table[i][j] <<" ";
        std::cout << "\n";
    }
}

// compute the average value of the intersection between
// the list of neighbours of two configurations
double Configuration::neighbour_overlap(Configuration b, bool sorting){
    double sum=0;
    if (sorting==false)
    {
        for (unsigned int i = 0; i < this->num_particles; ++i)
        {
            vector <int> common;
            set_intersection(this->neighbour_table[i].begin(), this->neighbour_table[i].end(), b.neighbour_table[i].begin(), b.neighbour_table[i].end(), back_inserter(common));
            
            // std::cout<<"particle "<<i<<std::endl;
            // for (int p = 0; p < common.size(); ++p)
            // {
            //     std::cout<<common[p]<<" ";
            // }
            // std::cout<<std::endl;
            
            // std::cout<<"==>"<<common.size()<<" "<<neighbour_table
            sum+=common.size()/(double)neighbour_table[i].size();
        }
    }
    else
    {
        // in case the neighbours are not sorted...
        for (unsigned int i = 0; i < this->num_particles; ++i)
        {
            vector <int> common;
            sort(this->neighbour_table[i].begin(), this->neighbour_table[i].end());
            sort(b.neighbour_table[i].begin(), b.neighbour_table[i].end());
            
            set_intersection(this->neighbour_table[i].begin(), this->neighbour_table[i].end(), b.neighbour_table[i].begin(), b.neighbour_table[i].end(), back_inserter(common));
            
            sum+=common.size()/neighbour_table[i].size();
        }
    }
    
    return sum/this->num_particles;
}

const vector<double>& Configuration::radial_distribution(unsigned int num_bins, double bin_width)
{
    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);
    if (this->g.size() == num_bins && this->g_bin_width == bin_width) return this->g;
    
    constexpr unsigned int d = 3;
    const double r_max = num_bins*bin_width;
    const double r_max_squ = pow(r_max, 2.);
    
    // Declare these to make finding particles code more legible.
    const ParticleIndex* id;
    const double* r1;
    const double* r2;
    // Calculation variables defined here for efficiency.
    double delta;
    double delta_r_squ;
    vector<unsigned int> bin_count(num_bins); // should automatically initialise to zeros.
    unsigned int bin;
    for (unsigned int i=0; i < (this->num_particles-1); ++i)
    {
        id = &this->particle_table[i];
        r1 = &this->particles[id->species][id->index];
        for (unsigned int j=i+1; j < this->num_particles; ++j)
        {
            id = &this->particle_table[j];
            r2 = &this->particles[id->species][id->index];
            
            delta_r_squ = 0.0;
            for (unsigned int c = 0; c < d; ++c)
            {
                delta = this->apply_boundaries(r1[c]-r2[c], c);
                delta_r_squ += delta*delta;
            }
            //delta_r = sqrt(delta_r);
            //this->compute_differences(deltas);
            //this->periodic_boundaries(deltas);
            
            // double r=norm(deltas); // calc distance
            
            if (delta_r_squ < r_max_squ) //within assigned distance from reference
            {
                bin = static_cast<unsigned int> ( sqrt(delta_r_squ)/bin_width );
                bin_count[bin] += 2;
            }
        }
    }
    
    this->g = vector<double>(num_bins);

    double number_density=this->num_particles/this->get_volume();

    for ( unsigned int i=0; i < num_bins; ++i )
    {
        this->g[i] = static_cast<double> ( bin_count[i] );
        //normalise
        double vol=((i+1)*(i+1)*(i+1)-i*i*i)*bin_width*bin_width*bin_width; //needs to be in 3D
        double nid=(4./3.)*M_PI*vol*number_density;//*this->number_density; //3D
        this->g[i]/=nid*this->num_particles; //scale         
    }
    
    return this->g;
}

