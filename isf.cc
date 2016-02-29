
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

