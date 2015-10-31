#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
using namespace std;

#include "configuration.h"
#include "utilities.h"


Configuration::Configuration() : num_particles(0), particles(0)
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
    constexpr int d = 3;

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
        for (unsigned int c = 0; c < d; ++c) in >> x[c];
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
            if (index_list.size() > this->dispersity.size())
                throw Exception(__PRETTY_FUNCTION__, ": more species found than expected",
                                ": found (at least)=", index_list.size(), ", expected=", this->dispersity.size());
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

        // Bookkeeping so every particle has a unique id (based on their order in the xyz file).
        this->particle_table[n].species = species_index;
        this->particle_table[n].index = count[species_index]-1;
    }
}

void Configuration::read_atom(string path)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_atom(in);
    in.close();
}

void Configuration::read_atom(istream& in)
{
    constexpr unsigned int d = 3;

    if (this->particles.size()) throw Exception(__PRETTY_FUNCTION__, ": attempting to load new configuration into a preexisting configuration");

    string line;
    double tmp;

    // The first two lines give the time (including a header line).
    getline(in, line);
    getline(in, line);
    //int t = stoi(line);

    // The next two give the number of atoms (including a header line).
    getline(in, line);
    getline(in, line);
    this->num_particles = stoi(line);
    this->particle_table = vector<ParticleIndex>(this->num_particles);

    // The next lines give the domain size (including a header line).
    getline(in, line);
    this->boundaries = vector<double>(d);
    for (unsigned int c = 0; c < d; ++c)
    {
        in >> tmp >> tmp;
        this->boundaries[c] = tmp;
    }
    getline(in, line); // finish loading the rest of this line.

    // Check the units for the particle positons: they may be stored in scaled units in which case we have to rescale the coordinates later on.
    getline(in, line);
    bool scaled_coordinates = line.find("s");

    // Declare all variables outside of the loop for optimisation.
    // NB: this could be unnecessary as this is not meant to be a fast function.
    int species, index;
    int species_index;
    map<int, int> index_list;
    map<int, int>::iterator found;
    vector< vector<double> > positions;
    double x[d];
    // Read particle data.
    for (unsigned int n = 0; n < this->num_particles; ++n)
    {
        in >> index >> species;
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
        for (unsigned int c = 0; c < d; ++c) in >> x[c];
        if (scaled_coordinates) for (unsigned int c = 0; c < d; ++c) x[c] *= this->boundaries[c];
        positions[species_index].insert(positions[species_index].end(), x, x+d);

        // Bookkeeping so every particle has a unique id.
        this->particle_table[index-1].species = species_index;
        this->particle_table[index-1].index = (positions[species_index].size()/d)-1;
    }

    // Finish off the current line (so we start at the beginning of the next configuration if there is more than one).
    getline(in, line);

    // Save the coordinate data as static-length species, using the minimum storage space.
    for (auto it = index_list.begin(); it != index_list.end(); ++it)
    {
        species_index = it->second;
        positions[species_index].resize( positions[species_index].size() );
        this->particles.push_back( Species<d>(positions[species_index]) );
        this->dispersity.push_back( positions[species_index].size()/d );
    }
}

void Configuration::read_atom(string path, const Configuration& ref_config)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_atom(in, ref_config);
    in.close();
}

void Configuration::read_atom(istream& in, const Configuration& ref_config)
{
    constexpr unsigned int d = 3;

    if (this->particles.size()) throw Exception(__PRETTY_FUNCTION__, ": attempting to load new configuration into a preexisting configuration");

    string line;
    double tmp;

    // The first two lines give the time (including a header line).
    getline(in, line);
    getline(in, line);
    //int t = stoi(line);

    // The next two give the number of atoms (including a header line).
    getline(in, line);
    getline(in, line);
    this->num_particles = stoi(line);

    if (this->num_particles != ref_config.num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to load configuration (N=", this->num_particles, ") which does not match reference configuration (N=", ref_config.num_particles, ")");

    this->particle_table = vector<ParticleIndex>(this->num_particles);

    // The next lines give the domain size (including a header line).
    getline(in, line);
    this->boundaries = vector<double>(d);

    for (unsigned int c = 0; c < d; ++c)
    {
        in >> tmp >> tmp;
        this->boundaries[c] = tmp;
        if (this->boundaries[c] != ref_config.boundaries[c]) throw Exception(__PRETTY_FUNCTION__, ": attempting to load configuration whose boundaries do not match reference configuration");
    }

    getline(in, line); // finish loading the rest of this line.

    // Preallocate data types for the different species.
    this->dispersity = vector<unsigned int>(ref_config.dispersity);
    for (auto it = this->dispersity.begin(); it != this->dispersity.end(); ++it)
        this->particles.push_back( Species<d>(*it) );

    // Check the units for the particle positons: they may be stored in scaled units in which case we have to rescale our coordinates later on.
    getline(in, line);
    bool scaled_coordinates = line.find("s");

    // Declare all variables outside of the loop for optimisation.
    unsigned int species, index;
    vector<unsigned int> count( this->dispersity.size() );
    // Read the particle positions.
    for (unsigned int n = 0; n < this->num_particles; ++n)
    {
        in >> index >> species;
        if (species > this->dispersity.size())
            throw Exception(__PRETTY_FUNCTION__, ": more species found than expected",
                            ": found (at least)=", species, ", expected=", this->dispersity.size());
        species--;

        // Check this new particle does not violate the prescribed dispersity.
        count[species] += 1;
        if (count[species] > this->dispersity[species])
            throw Exception(__PRETTY_FUNCTION__, ": more particles found of species ", species, " than expected",
                            ": found (at least)=", count[species], ", expected=", this->dispersity[species]);

        // Get the coordinates.
        for (unsigned int c = 0; c < d; c++) in >> this->particles[species]( count[species]-1, c );
        if (scaled_coordinates)
        {
            for (unsigned int c = 0; c < d; ++c) this->particles[species]( count[species]-1, c ) *= this->boundaries[c];
        }

        // Bookkeeping so every particle has a unique id.
        this->particle_table[index-1].species = species;
        this->particle_table[index-1].index = count[species]-1;
    }

    // Finish off the current line (so we start at the beginning of the next configuration if there is more than one).
    getline(in, line);
}


void Configuration::read_neighbours(std::string path)
{
    std::ifstream file(path);
    if (!file.good()) std::cerr << "ERROR: the file " << path << " does not exist\n Forced exit.\n";
    std::vector< std::vector<int> >table;
    for (std::string line; std::getline(file, line); )
    {
        std::stringstream stream;
        stream.str(line);

        int count=0;
        std::vector<int> neighs;
        for (std::string s; std::getline(stream, s,' '); ){
            // skip the first element
            if (count > 0)
            {
                neighs.push_back(StringToNum<int>(s));
            }
            count++;
        }
        // cout<<neighs.size()<<"\n";
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
        out << "species " << i << ": " << sp.size() << " particles\n";
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
double Configuration::neighbour_overlap(Configuration& b, bool sorting)
{
    double sum=0;
    if (sorting==false)
    {
        for (unsigned int i = 0; i < this->num_particles; ++i)
        {
            vector <int> common;
            set_intersection(this->neighbour_table[i].begin(), this->neighbour_table[i].end(), b.neighbour_table[i].begin(), b.neighbour_table[i].end(), back_inserter(common));

            // cout << "particle " << i << "\n";
            // for (int p = 0; p < common.size(); ++p)
            // {
            //     cout << common[p] << " ";
            // }
            // cout << "\n";

            // cout << "==>" << common.size() << " " << neighbour_table;
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

vector<double> Configuration::radial_distribution(unsigned int num_bins, double bin_width) const
{
    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);

    vector<double> g(num_bins);
    this->cumulative_radial_distribution(g, bin_width);
    return g;
}

vector<double> Configuration::radial_distribution(unsigned int species, unsigned int num_bins, double bin_width) const
{
    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);
    if (species >= this->dispersity.size()) throw Exception(__PRETTY_FUNCTION__, ": invalid first species=", species);

    vector<double> g(num_bins);
    this->cumulative_radial_distribution(species, g, bin_width);
    for (unsigned int i = 0; i < num_bins; ++i)
        g[i] = g[i] * (1.0*this->num_particles/this->dispersity[species])*(1.0*this->num_particles/this->dispersity[species]);

    return g;
}

vector<double> Configuration::radial_distribution(unsigned int species_a, unsigned int species_b, unsigned int num_bins, double bin_width) const
{
    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);
    if (species_a >= this->dispersity.size()) throw Exception(__PRETTY_FUNCTION__, ": invalid first species=", species_a);
    if (species_b >= this->dispersity.size()) throw Exception(__PRETTY_FUNCTION__, ": invalid second species=", species_b);

    vector<double> g(num_bins);
    this->cumulative_radial_distribution(species_a, species_b, g, bin_width);
    return g;
}

void Configuration::cumulative_radial_distribution(vector<double>& g_total, double bin_width) const
{
    const unsigned int num_bins = g_total.size();

    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);

    for (unsigned int species_a = 0; species_a < this->dispersity.size(); ++species_a)
        for (unsigned int species_b = species_a; species_b < this->dispersity.size(); ++species_b)
            this->cumulative_radial_distribution(species_a, species_b, g_total, bin_width);
}

void Configuration::cumulative_radial_distribution(unsigned int species, vector<double>& g_total, double bin_width) const
{
    const unsigned int num_bins = g_total.size();

    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);
    if (species >= this->dispersity.size()) throw Exception(__PRETTY_FUNCTION__, ": invalid species=", species);

    this->cumulative_radial_distribution(species, species, g_total, bin_width);
}

void Configuration::cumulative_radial_distribution(unsigned int species_a, unsigned int species_b, vector<double>& g_total, double bin_width) const
{
    constexpr unsigned int d = 3;
    const unsigned int num_bins = g_total.size();
    const double r_max = num_bins*bin_width;
    const double r_max_squ = r_max*r_max;

    if (!this->num_particles) throw Exception(__PRETTY_FUNCTION__, ": attempting to compute g(r) on an empty configuration");
    if (!num_bins) throw Exception(__PRETTY_FUNCTION__, ": invalid num_bins=", num_bins);
    //if (!this->bin_width > 0.) throw Exception(__PRETTY_FUNCTION__, ": invalid bin_width=", bin_width);
    if (species_a >= this->dispersity.size()) throw Exception(__PRETTY_FUNCTION__, ": invalid species=", species_a);
    if (species_b >= this->dispersity.size()) throw Exception(__PRETTY_FUNCTION__, ": invalid species=", species_b);

    // Declare these to make the 'find-particle' code more legible.
    const double* r1;
    const double* r2;
    // Calculation variables defined here for efficiency.
    double delta;
    double delta_r_squ;
    vector<unsigned int> bin_count(num_bins); // should automatically initialise to zeros.
    unsigned int bin;

    // We have to iterate slightly differently if comparing the same species to avoid double counting.
    const unsigned int i_start = 0;
    const unsigned int i_end = this->dispersity[species_a] - (species_a == species_b);
    const unsigned int j_end = this->dispersity[species_b];

    for (unsigned int i = i_start; i < i_end; ++i)
    {
        r1 = &this->particles[species_a][i];

        const unsigned int j_start = (species_a == species_b) ? i+1 : 0;
        for (unsigned int j = j_start; j < j_end; ++j)
        {
            r2 = &this->particles[species_b][j];

            delta_r_squ = 0.0;
            for (unsigned int c = 0; c < d; ++c)
            {
                delta = this->apply_boundaries(r1[c]-r2[c], c);
                delta_r_squ += delta*delta;
            }

            if (delta_r_squ < r_max_squ) //within assigned distance from reference
            {
                bin = static_cast<unsigned int> ( sqrt(delta_r_squ)/bin_width );
                bin_count[bin] += 2; // 2 particles contribute
            }
        }
    }

    // We need the global number density for normalisation.
    // If computing the density for a single species then further renormalisation using the species number density may be needed.
    const double number_density = this->num_particles/this->get_volume();

    for ( unsigned int i = 0; i < num_bins; ++i )
    {
        // This normalisation is only valid in 3d:
        double vol=((i+1)*(i+1)*(i+1)-i*i*i)*bin_width*bin_width*bin_width;
        double nid=(4./3.)*M_PI*vol*number_density;
        g_total[i] += bin_count[i]/(nid*this->num_particles);
    }
}

vector<double> Configuration::msd_isf(const Configuration& b, const double q) const
{
    constexpr unsigned int d = 3;
    vector<double> msd_isf(d+1);
    this->cumulative_msd_isf(msd_isf, b, q);
    return msd_isf;
}

void Configuration::cumulative_msd_isf(vector<double>& msd_isf_total, const Configuration& b, const double q) const
{
    constexpr unsigned int d = 3;

    const ParticleIndex* lookup;
    const double* ra;
    const double* rb;
    double delta;

    double dr;
    vector<double> dr_squ(d);

    for (unsigned int i = 0; i < this->num_particles; ++i)
    {
        lookup = &this->particle_table[i];
        ra = &this->particles[lookup->species][lookup->index];

        lookup = &b.particle_table[i];
        rb = &b.particles[lookup->species][lookup->index];

        dr = 0.;
        for (unsigned int c = 0; c < d; ++c)
        {
            delta = this->apply_boundaries(rb[c]-ra[c], c);
            dr_squ[c] = delta*delta;
            dr += dr_squ[c];
            msd_isf_total[c] += dr_squ[c];
        }
        dr = sqrt(dr);
        msd_isf_total[d] += sin(q*dr)/(q*dr);
    }
}
