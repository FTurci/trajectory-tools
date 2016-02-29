#include "configuration.h"
#include "utilities.h"

#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;


/*** The abstract base class for all configuration types. ***/


BaseConfiguration::BaseConfiguration(const Container& copy)
    : Container(copy)
{
}

void BaseConfiguration::write_coords(string path) const
{
    ofstream out(path);
    if (!out) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->write_coords(out);
    out.close();
}

void BaseConfiguration::write_coords(ostream& out) const
{
    const ParticleIndex* lookup;
    const double* r;

    for (unsigned int i = 0; i < this->diversity.system_size; ++i)
    {
        lookup = &this->particle_table[i];
        r = this->position(lookup->species, lookup->index);
        for (unsigned int c = 0; c < d; ++c)
            out << "  " << this->apply_boundaries(r[c]-this->origin[c], c);
        out << "\n";
    }
}

void BaseConfiguration::write_xyz(string path) const
{
    ofstream out(path);
    if (!out) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->write_xyz(out);
    out.close();
}

void BaseConfiguration::write_xyz(ostream& out) const
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
        r = this->position(lookup->species, lookup->index);
        for (unsigned int c = 0; c < d; ++c)
            out << "  " << this->apply_boundaries(r[c]-this->origin[c], c);
        out << "\n";
    }
}

array<double,Container::d> BaseConfiguration::pseudo_box() const
{
    array<double,d> box;
    const ParticleIndex* lookup;
    const double* r;

    for (unsigned int i = 0; i < this->size(); ++i)
    {
        lookup = &this->particle_table[i];
        r = this->position(lookup->species, lookup->index);
        for (unsigned int c = 0; c < d; ++c) box[c] = max(box[c], r[c]);
    }

    return box;
}

void BaseConfiguration::read_xyz(string path)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_xyz(in);
    in.close();
}

void BaseConfiguration::read_xyz(istream& in)
{
    if (this->is_empty())
        throw Exception(__PRETTY_FUNCTION__,
                        ": attempting to load configuration into an unknown data structure");

    // The first line states the number of particles.
    string line;
    getline(in, line);
    if (stoul(line) != this->size())
        throw Exception(__PRETTY_FUNCTION__, ": attempting to read configuration with ", stoi(line), " particles when configuration expects ", this->diversity.system_size);

    // The comment line can be ignored.
    getline(in, line);

    // Read the coordinates.
    string species;
    int species_index;
    map<string, unsigned int>::iterator found;
    vector<unsigned int> count( this->num_species() );
    for (unsigned int n = 0; n < this->size(); ++n)
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
        {
            double temp;
            in >> temp;
            this->position(species_index, count[species_index]-1)[c] = temp + this->origin[c];
        }

        // Bookkeeping so every particle has a unique id (based on their order in the xyz file).
        this->particle_table[n].species = species_index;
        this->particle_table[n].index = count[species_index]-1;
    }

    // // Check we found the right number of particles.
    // for (unsigned int index = 0; index < this->num_species(); ++index)
    // {
    //     int expected = this->diversity.species[index].second;
    //     int actual = count[index];
    //     if (actual != expected)
    //         throw Exception(__PRETTY_FUNCTION__, ": not enough particles read of species ",
    //                         this->diversity.species[index].first,
    //                         ", expected=", expected, " but found=", actual);
    // }
}


/*** View class for examining subconfigurations. ***/


ConfigurationView::ConfigurationView(const BaseConfiguration* source)
    : BaseConfiguration( *dynamic_cast<const Container*>(source) )
{
    this->diversity = Diversity(source->summary());
    this->particle_table.reserve(this->size());
    this->positions = std::vector< std::vector<double*> >(this->num_species());
    for (unsigned int i = 0; i < this->num_species(); ++i)
    {
        this->positions[i].reserve(this->diversity.species[i].second);
        this->diversity.species[i].second = 0;
    }
}


/*** Main configuration class, containing all the entire system ***/


Configuration::Configuration() : particles(0)
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

void Configuration::read_xyz(string path, const BaseConfiguration* ref_config)
{
    this->read_xyz(path, ref_config->summary());
}

void Configuration::read_xyz(string path, const Diversity& ref_div)
{
    ifstream in(path);
    if (!in) throw Exception(__PRETTY_FUNCTION__, ": could not open ", path);
    this->read_xyz(in, ref_div);
    in.close();
}

void Configuration::read_xyz(istream& in, const BaseConfiguration* ref_config)
{
    this->read_xyz(in, ref_config->summary());
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
    for (unsigned int n = 0; n < this->size(); ++n)
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

void Configuration::read_xyz(istream& in, const Diversity& ref_div)
{
    if (!this->is_empty())
        throw Exception(__PRETTY_FUNCTION__,
                        ": attempting to load new configuration into a preexisting configuration");
    this->diversity = Diversity(ref_div);
    const unsigned int num_species = this->num_species();

    // The first line states the number of particles.
    string line;
    getline(in, line);
    if (stoul(line) != this->size())
        throw Exception(__PRETTY_FUNCTION__, ": attempting to read configuration with ", stoi(line), " particles when configuration expects ", this->diversity.system_size);

    // The comment line is not standardised. It may contain the box size however.
    getline(in, line);
    bool found_box = false;
    {
        stringstream ss(line);
        string first_word;
        ss >> first_word;
        if (first_word == "box")
        {
            for (unsigned int c = 0; c < d; ++c)
                ss >> this->boundaries[c];
            found_box = true;
        }
    }

    // Preallocate the data types.
    this->particles.reserve(this->num_species());
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


/*Configuration Configuration::subsystem(double distance) const
{
    double center[d] = {0.0, 0.0, 0.0};  /// needs to be generalised to d-dimensions
    return this->subsystem(distance, center);
}

Configuration Configuration::subsystem(double distance, unsigned int central_particle) const
{
    auto lookup = &this->particle_table[central_particle];
    const double* center = &this->particles[lookup->species][lookup->index];
    return this->subsystem(distance, center);
    }*/


/*Configuration Configuration::exclude_subsystem(double distance) const
{
    double center[d] = {0.0, 0.0, 0.0};  /// needs to be generalised to d-dimensions
    return this->exclude_subsystem(distance, center);
}

Configuration Configuration::exclude_subsystem(double distance, unsigned int central_particle) const
{
    auto lookup = &this->particle_table[central_particle];
    const double* center = &this->particles[lookup->species][lookup->index];
    return this->exclude_subsystem(distance, center);
    }*/


/*void Configuration::set_origin(array<double,d> new_origin)
{
    for (unsigned int i = 0; i < this->diversity.system_size; ++i)
    {
        auto lookup = &this->particle_table[i];
        double* r = this->position(lookup->species, lookup->index);
        for (unsigned int c = 0; c < d; ++c)
            r[c] = apply_boundaries(r[c] - new_origin[c], c);
    }

    }*/
