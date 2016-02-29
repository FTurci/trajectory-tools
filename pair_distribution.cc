
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
