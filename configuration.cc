#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

#include "configuration.h"
#include "utilities.h"


Configuration::Configuration() : num_particles(0), particles(0)
{
}


/// XYZ files.


void Configuration::read_xyz(string path)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_xyz(in);
    in.close();
}

void Configuration::read_xyz(string path, const Configuration& ref_config)
{
    this->read_xyz(path, ref_config.diversity);
}

void Configuration::read_xyz(string path, const Diversity& ref_div)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_xyz(in, ref_div);
    in.close();
}

void Configuration::read_xyz(istream& in)
{
    if (this->particles.size()) throw Exception(__PRETTY_FUNCTION__, ": attempting to load new configuration into a preexisting configuration");

    // The first line states the number of particles.
    string line;
    getline(in, line);
    this->diversity.system_size = stoi(line);
    this->particle_table = vector<ParticleIndex>(this->diversity.system_size);

    // The comment line is not standardised. It may contain the box size however.
    getline(in, line);
    bool found_box = false;
    {
        stringstream ss(line);
        string first_word;
        ss >> first_word;
        if (first_word == "box")
        {
            this->boundaries = vector<double>(d);
            for (unsigned int c = 0; c < d; ++c)
                ss >> this->boundaries[c];
            found_box = true;
        }
    }

    // Read the coordinates.
    string species;
    int species_index;
    map<string, unsigned int>::iterator found;
    vector< vector<double> > positions;
    double x[d];
    for (unsigned int n = 0; n < this->diversity.system_size; ++n)
    {
        in >> species;
        // If this is a new species then we have to add new data structures for it.
        found = this->diversity.species_map.find(species);
        if (found == this->diversity.species_map.end())
        {
            positions.push_back( vector<double>() );
            species_index = this->diversity.species_map.size();
            this->diversity.species_map[species] = species_index;
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
    const unsigned int num_species = positions.size();
    this->particles.reserve(num_species);
    for (auto it = positions.begin(); it != positions.end(); ++it)
        this->particles.push_back(*it);
    // And the summary information.
    this->diversity.species = vector< pair<string,unsigned int> >(num_species);
    for (auto it = this->diversity.species_map.begin(); it != this->diversity.species_map.end(); ++it)
    {
        species_index = it->second;
        this->diversity.species[species_index] = pair<string,unsigned int> (it->first,positions[species_index].size()/d);
    }

    // Xyz files do not store the box size by default so we may have to approximate it.
    if (!found_box) this->boundaries = this->pseudo_box();;
}

void Configuration::read_xyz(std::istream& in, const Diversity& ref_div)
{
    if (this->particles.size())
        throw Exception(__PRETTY_FUNCTION__,
                        ": attempting to load new configuration into a preexisting configuration");
    this->diversity = Diversity(ref_div);
    const unsigned int num_species = this->diversity.species.size();

    // The first line states the number of particles.
    string line;
    getline(in, line);
    if (stoul(line) != this->diversity.system_size)
        throw Exception(__PRETTY_FUNCTION__, ": attempting to read configuration with ", stoi(line),
                        " particles when configuration expects ", this->diversity.system_size);

    // The comment line is not standardised. It may contain the box size however.
    getline(in, line);
    bool found_box = false;
    {
        stringstream ss(line);
        string first_word;
        ss >> first_word;
        if (first_word == "box")
        {
            this->boundaries = vector<double>(d);
            for (unsigned int c = 0; c < d; ++c)
                ss >> this->boundaries[c];
            found_box = true;
        }
    }

    // Preallocate the data types.
    this->particles.reserve(num_species);
    for (auto it = this->diversity.species.begin(); it != this->diversity.species.end(); ++it)
        this->particles.push_back( Species<d>(it->second) );
    this->particle_table = vector<ParticleIndex>(this->diversity.system_size);

    // Read the coordinates.
    string species;
    int species_index;
    map<string, unsigned int>::iterator found;
    vector<unsigned int> count( num_species );
    for (unsigned int n = 0; n < this->diversity.system_size; ++n)
    {
        in >> species;

        found = this->diversity.species_map.find(species);
        if (found == this->diversity.species_map.end())
            throw Exception(__PRETTY_FUNCTION__, ": unrecognised species ", species);
        else species_index = found->second;

        // Check this new particle does not violate the prescribed dispersity.
        count[species_index] += 1;
        if (count[species_index] > this->diversity.species[species_index].second)
            throw Exception(__PRETTY_FUNCTION__, ": more particles found of species ", species, " than expected: found (at least)=", count[species_index], ", expected=", this->diversity.species[species_index].second);

        // Get the coordinates.
        for (unsigned int c = 0; c < d; ++c)
            in >> this->particles[species_index]( count[species_index]-1, c );

        // Bookkeeping so every particle has a unique id (based on their order in the xyz file).
        this->particle_table[n].species = species_index;
        this->particle_table[n].index = count[species_index]-1;
    }

    // Xyz files do not store the box size by default so we may have to approximate it.
    if (!found_box) this->boundaries = this->pseudo_box();;
}

void Configuration::write_xyz(string path) const
{
    ofstream out(path);
    if (!out) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->write_xyz(out);
    out.close();
}

void Configuration::write_xyz(ostream& out) const
{
    out << this->diversity.system_size << "\n";
    // The next line is a comment, so it does not strictly matter.
    // We may as well put the box size here:
    out << "box";
    for (unsigned int c = 0; c < d; ++c) out << " " << this->boundaries[c];
    out << "\n";

    const ParticleIndex* lookup;
    const double* r;

    for (unsigned int i = 0; i < this->diversity.system_size; ++i)
    {
        lookup = &this->particle_table[i];
        out << this->diversity.species[lookup->species].first;
        r = &this->particles[lookup->species][lookup->index];
        for (unsigned int c = 0; c < d; ++c) out << "  " << r[c];
        out << "\n";
    }
}


/// LAMMPS files



void Configuration::read_atom(string path)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_atom(in);
    in.close();
}

void Configuration::read_atom(istream& in)
{
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
    for (unsigned int i = 0; i < this->diversity.species.size(); ++i)
    {
        const Species<d>& sp = this->particles[i];
        out << "species " << i+1 << " (" << this->diversity.species[i].first << "): ";
        out << sp.size() << " particles\n";
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


vector<double> Configuration::pseudo_box() const
{
    vector<double> box(d, 0.0);
    const ParticleIndex* lookup;
    const double* r;

    for (unsigned int i = 0; i < this->num_particles; ++i)
    {
        lookup = &this->particle_table[i];
        r = &this->particles[lookup->species][lookup->index];
        for (unsigned int c = 0; c < d; ++c) box[c] = max(box[c], r[c]);
    }

    return box;
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
    vector<double> msd_isf(d+1);
    this->cumulative_msd_isf(msd_isf, b, q);
    return msd_isf;
}

void Configuration::cumulative_msd_isf(vector<double>& msd_isf_total, const Configuration& b, const double q) const
{
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


Configuration Configuration::subsystem(double distance) const
{
    double center[d] = {0.0, 0.0, 0.0};  /// needs to be generalised to d-dimensions
    return this->subsystem(distance, center);
}

Configuration Configuration::subsystem(double distance, unsigned int central_particle) const
{
    auto lookup = &this->particle_table[central_particle];
    const double* center = &this->particles[lookup->species][lookup->index];
    return this->subsystem(distance, center);
}

Configuration Configuration::subsystem(double distance, const double* from) const
{
    Configuration sub;
    sub.boundaries = this->boundaries;
    sub.diversity = Diversity(this->diversity);
    sub.diversity.system_size = 0;
    for (auto it = sub.diversity.species.begin(); it != sub.diversity.species.end(); ++it)
        it->second = 0;

    vector< vector<double> > positions(this->diversity.species.size());
    sub.particle_table.reserve(this->diversity.system_size);

    double delta, dr_squ;
    const double distance_squ = distance*distance;

    for (unsigned int i = 0; i < this->diversity.system_size; ++i)
    {
        auto lookup = &this->particle_table[i];
        const double* r = &this->particles[lookup->species][lookup->index];

        dr_squ = 0.;
        for (unsigned int c = 0; c < d; ++c)
        {
            delta = this->apply_boundaries(r[c]-from[c], c);
            dr_squ += delta*delta;
        }

        if (dr_squ < distance_squ)
        {
            positions[lookup->species].insert(positions[lookup->species].end(), r, r+d);

            // Update the bookkeeping.
            sub.particle_table.push_back(ParticleIndex());
            sub.particle_table[sub.diversity.system_size].species = lookup->species;
            sub.particle_table[sub.diversity.system_size].index = sub.diversity.species[lookup->species].second;
            sub.diversity.species[lookup->species].second++;
            sub.diversity.system_size++;
        }
    }

    sub.particle_table.shrink_to_fit();
    for (auto it = positions.begin(); it != positions.end(); ++it)
        sub.particles.push_back(*it);

    return sub;
}

Configuration Configuration::exclude_subsystem(double distance) const
{
    double center[d] = {0.0, 0.0, 0.0};  /// needs to be generalised to d-dimensions
    return this->exclude_subsystem(distance, center);
}

Configuration Configuration::exclude_subsystem(double distance, unsigned int central_particle) const
{
    auto lookup = &this->particle_table[central_particle];
    const double* center = &this->particles[lookup->species][lookup->index];
    return this->exclude_subsystem(distance, center);
}

Configuration Configuration::exclude_subsystem(double distance, const double* from) const
{
    Configuration sub;
    sub.boundaries = this->boundaries;
    sub.diversity = Diversity(this->diversity);
    sub.diversity.system_size = 0;
    for (auto it = sub.diversity.species.begin(); it != sub.diversity.species.end(); ++it)
        it->second = 0;

    vector< vector<double> > positions(this->diversity.species.size());
    sub.particle_table.reserve(this->diversity.system_size);

    double delta, dr_squ;
    const double distance_squ = distance*distance;

    for (unsigned int i = 0; i < this->diversity.system_size; ++i)
    {
        auto lookup = &this->particle_table[i];
        const double* r = &this->particles[lookup->species][lookup->index];

        dr_squ = 0.;
        for (unsigned int c = 0; c < d; ++c)
        {
            delta = this->apply_boundaries(r[c]-from[c], c);
            dr_squ += delta*delta;
        }

        if (!(dr_squ < distance_squ))
        {
            positions[lookup->species].insert(positions[lookup->species].end(), r, r+d);

            // Update the bookkeeping.
            sub.particle_table.push_back(ParticleIndex());
            sub.particle_table[sub.diversity.system_size].species = lookup->species;
            sub.particle_table[sub.diversity.system_size].index = sub.diversity.species[lookup->species].second;
            sub.diversity.species[lookup->species].second++;
            sub.diversity.system_size++;
        }
    }

    sub.particle_table.shrink_to_fit();
    for (auto it = positions.begin(); it != positions.end(); ++it)
        sub.particles.push_back(*it);

    return sub;
}

void Configuration::set_origin(array<double,d> new_origin)
{
    for (unsigned int i = 0; i < this->diversity.system_size; ++i)
    {
        auto lookup = &this->particle_table[i];
        double* r = &this->particles[lookup->species][lookup->index];
        for (unsigned int c = 0; c < d; ++c)
            r[c] = apply_boundaries(r[c] - new_origin[c], c);
    }

}
