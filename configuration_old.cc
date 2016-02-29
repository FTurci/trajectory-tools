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
