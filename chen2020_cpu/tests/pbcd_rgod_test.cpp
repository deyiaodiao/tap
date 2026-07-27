#include "PBCD_Algorithm.h"

#include <cassert>
#include <cmath>

int main()
{
    const double gap = TAP_PBCD::ComputeODRelativeGap(100.0, 10.0, 1080.0);
    assert(std::abs(gap - (1.0 - 1000.0 / 1080.0)) < 1e-15);
    assert(TAP_PBCD::ComputeODRelativeGap(100.0, 10.0, 1000.0) == 0.0);
    assert(TAP_PBCD::ComputeODRelativeGap(100.0, 10.0, 0.0) == 0.0);
    return 0;
}
