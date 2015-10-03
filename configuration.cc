#include <iostream>
#include <fstream>
#include <sstream>  
#include <algorithm>
#include <map>
using namespace std;

#include "configuration.h"
#include "utilities.h"


Configuration::Configuration() : numParticles(0)
{

}

void Configuration::read_xyz(string path)
{
    const int d = 3;
    
    if (this->particles.size()) throw Exception(__PRETTY_FUNCTION__, ": attempting to load new configuration into a preexisting configuration");
    
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    
    // The first line states the number of particles.
    string line;
    getline(in, line);
    this->numParticles = stoi(line);
    
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
    for (unsigned int n = 0; n < this->numParticles; ++n)
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
    }
    
    in.close();
    
    // Save the coordinate data as static-length species, using the minimum storage space.
    for (auto it = index_list.begin(); it != index_list.end(); ++it)
    {
        species_index = it->second;
        positions[species_index].resize( positions[species_index].size() );
        this->particles.push_back( Species<d>(positions[species_index]) );
        this->dispersity.push_back( positions[species_index].size()/d );
    }
    
    /*********** Debug: print out every species' coordinate data. *************/
    for (unsigned int i = 0; i < this->particles.size(); ++i)
    {
        const Species<3>& sp = this->particles[i];
        cout << "species " << i << ": " << sp.size() << " particles" << endl;
        for (unsigned int n = 0; n < sp.size(); ++n)
        {
            for (unsigned int c = 0; c < d; ++c)
                cout << "  " << sp(n,c);
            cout << endl;
        }
    }
    /**********                   <\Debug>                ********************/
}

void Configuration::read_xyz(string path, const vector<unsigned int>& species_distribution)
{
    const int d = 3;
    
    if (this->particles.size())
        throw Exception(__PRETTY_FUNCTION__,
                        ": attempting to load new configuration into a preexisting configuration");
    
    // Make sure the given distribution is compatible with any other data structures we have.
    unsigned int n = 0;
    for (auto it = species_distribution.begin(); it != species_distribution.end(); ++it) n += *it;
    if (!this->numParticles) this->numParticles = n;
    else
    {
        if (this->numParticles != n)
            throw Exception(__PRETTY_FUNCTION__, ": attempting to load incompatible configuration: ",
                            "new size=", n, " but expected size=", this->numParticles);
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
    
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    
    // The first line states the number of particles.
    string line;
    getline(in, line);
    if (stoul(line) != this->numParticles)
        throw Exception(__PRETTY_FUNCTION__, ": file ", path, " contains ", stoi(line),
                        " particles when configuration expects ", this->numParticles);
    
    // The comment line can just be discarded (for now, this may change in later revisions).
    getline(in, line);
    
    // Preallocate the data types.
    for (auto it = this->dispersity.begin(); it != this->dispersity.end(); ++it)
    {
        this->particles.push_back( Species<d>(*it) );
    }
    
    // Declare all variables outside of the loop for optimisation.
    string species;
    int species_index;
    map<string, int> index_list;
    map<string, int>::iterator found;
    vector<unsigned int> count( this->dispersity.size() );
    for (unsigned int n = 0; n < this->numParticles; ++n)
    {
        in >> species;
        // If this is a new species then we have to add new data structures for it.
        found = index_list.find(species);
        if (found == index_list.end())
        {
            species_index = index_list.size();
            index_list[species] = species_index;
            if (index_list.size() > species_distribution.size())
                throw Exception(__PRETTY_FUNCTION__, ": more species found than expected in ", path,
                                ": found (at least)=", index_list.size(), ", expected=", species_distribution.size());
        }
        else species_index = found->second;
        
        // Check this new particle does not violate the prescribed dispersity.
        count[species_index] += 1;
        if (count[species_index] > this->dispersity.size())
            throw Exception(__PRETTY_FUNCTION__, ": more particles found in of species ", species, " than expected in ",
                            path, ": found (at least)=", count[species_index], ", expected=", this->dispersity.size());
        
        // Get the coordinates.
        for (unsigned int c = 0; c < d; c++)
            in >> this->particles[species_index]( count[species_index]-1, c );
    }
    
    in.close();
    
    /*********** Debug: print out every species' coordinate data. *************/
    for (unsigned int i = 0; i < this->particles.size(); ++i)
    {
        const Species<d>& sp = this->particles[i];
        cout << "species " << i << ": " << sp.size() << " particles" << endl;
        for (unsigned int n = 0; n < sp.size(); ++n)
        {
            for (unsigned int c = 0; c < d; ++c)
                cout << "  " << sp(n,c);
            cout << endl;
        }
    }
    /**********                   <\Debug>                ********************/
}

void Configuration::read_neighbours(std::string filename)
{
    std::ifstream file(filename.c_str());
    if(!file.good()){std::cerr<<"ERROR: the file "<<filename<<" does not exist\n Forced exit.\n";}
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
    this->Npart=neighbour_table.size();

}

void Configuration::print_neighbours(int first, int last){
    for (int i = first; i < last; ++i)
    {
        for (unsigned int j = 0; j < neighbour_table[i].size(); ++j)
        {
            std::cout<<this->neighbour_table[i][j]<<" ";
        }
            std::cout<<std::endl;
    }

}
void Configuration::print_neighbours(){
    for (int i = 0; i < this->Npart; ++i)
    {
        for (unsigned int j = 0; j < neighbour_table[i].size(); ++j)
        {
            std::cout<<this->neighbour_table[i][j]<<" ";
        }
            std::cout<<std::endl;
    }

}

// compute the average value of the intersection between
// the list of neighbours of two configurations
double Configuration::neighbour_overlap(Configuration b, bool sorting){
    double sum=0;
    if(sorting==false){
        for (int i = 0; i < this->Npart; ++i)
        {   
            std::vector <int> common;
            std::set_intersection(this->neighbour_table[i].begin(), this->neighbour_table[i].end(), b.neighbour_table[i].begin(),  b.neighbour_table[i].end(), std::back_inserter(common));

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
    else{
        // in case the neighbours are not sorted...
        for (int i = 0; i < this->Npart; ++i)
        { 
            std::vector <int> common;
            std::sort(this->neighbour_table[i].begin(), this->neighbour_table[i].end());
            std::sort(b.neighbour_table[i].begin(), b.neighbour_table[i].end());

            std::set_intersection(this->neighbour_table[i].begin(), this->neighbour_table[i].end(), b.neighbour_table[i].begin(),  b.neighbour_table[i].end(), std::back_inserter(common));

            sum+=common.size()/neighbour_table[i].size();
        }
    }


    return sum/this->Npart;
}


void Configuration::radial_distr(int nbins,double biwidth){
    cout << nbins << endl;
    cout << biwidth << endl;
/*    double rmax=nbins*binwidth;
    int bin;
    this->g.resize(nbins);
    // dx,dy,dz...
    std::vector<double> deltas;
    for (int i=0; i<this->Npart-1; i++) 
    {
        for (int j=i+1; j<this->Npart; j++)
        {
            this->compute_differences(deltas);
            this->periodic_boundaries(deltas);

            double r=norm(deltas); // calc distance

            if (r<rmax) //within assigned distance from reference
                {   
                    bin= (int)(r/binwidth);
                    gr[bin]+=2;
                }
        }
    }

    for ( int i=0; i< nbins; ++i )
        {
            //normalise
            double vol=((i+1)*(i+1)*(i+1)-i*i*i)*binwidth*binwidth*binwidth; //needs to be in 3D
            double nid=(4./3.)*M_PI*vol*this->density; //3D
            g[i]/=nid*this->Npart; //scale         
        }
*/
}

