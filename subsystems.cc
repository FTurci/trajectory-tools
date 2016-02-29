#include "subsystems.h"

OutsideBubble::OutsideBubble(const Container& system,
              const array<double,3>& bubble_center,
              const double bubble_radius)
    : Container(system), radius_squ(bubble_radius*bubble_radius)
{
    for (int c = 0; c < d; ++c) center[c] = bubble_center[c] + this->origin[c];
}
