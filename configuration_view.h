#ifndef __CONFIGURATION_VIEW_H
#define __CONFIGURATION_VIEW_H

#include "configuration_base.h"


// Only supports canonical configurations for now, i.e. particle number N is fixed.
class ConfigurationView : public BaseConfiguration
{
public:
    friend class BaseConfiguration;
    ConfigurationView(const BaseConfiguration* source);

    template<class Criteria> ConfigurationView subsystem(Criteria criteria) const
    {
        return BaseConfiguration::subsystem<ConfigurationView>(criteria);
    }

protected:
    // Particle positions and their access.
    std::vector< std::vector<double*> > positions;
    inline const double* position(int s, int n) const
    {
        return this->positions[s][n];
    }
    inline double* position(int s, int n)
    {
        return this->positions[s][n];
    }
};

#endif
