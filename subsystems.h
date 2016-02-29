#ifndef __SUBSYSTEMS_H
#define __SUBSYSTEMS_H

#include "container.h"

// Here are a few useful functors for creating subsystem views for use with
// Configuration.subsystem(...) function.


class InsideBubble : public Container
{
public:
    inline InsideBubble(const Container& system,
                        const std::array<double,3>& bubble_center,
                        const double bubble_radius)
        : Container(system), radius_squ(bubble_radius*bubble_radius)
    {
        for (int c = 0; c < d; ++c) center[c] = bubble_center[c] + this->origin[c];
    }

    inline bool operator() (const double* r) const
    {
        double dr_squ = 0.;
        for (unsigned int c = 0; c < d; ++c)
        {
            double delta = this->apply_boundaries(r[c]-center[c], c);
            dr_squ += delta*delta;
        }
        return dr_squ < this->radius_squ;
    }

private:
    std::array<double,3> center;
    double radius_squ;
};

class OutsideBubble : public Container
{
public:
    inline OutsideBubble(const Container& system,
                  const std::array<double,3>& bubble_center,
                  const double bubble_radius)
        : Container(system), radius_squ(bubble_radius*bubble_radius)
    {
        for (int c = 0; c < d; ++c) center[c] = bubble_center[c] + this->origin[c];
    }

    inline bool operator() (const double* r) const
    {
        double dr_squ = 0.;
        for (unsigned int c = 0; c < d; ++c)
        {
            double delta = this->apply_boundaries(r[c]-center[c], c);
            dr_squ += delta*delta;
        }
        return !(dr_squ < this->radius_squ);
    }

private:
    std::array<double,3> center;
    double radius_squ;
};

#endif
