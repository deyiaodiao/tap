#include "PBCD_Algorithm.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <vector>

namespace
{
void CheckCoverage(std::size_t count, std::size_t blockSize)
{
    const auto blocks = TAP_PBCD::BuildConstantDistanceBlocks(count, blockSize);
    std::vector<std::size_t> flattened;
    for (const auto& block : blocks)
        flattened.insert(flattened.end(), block.begin(), block.end());
    std::sort(flattened.begin(), flattened.end());
    std::vector<std::size_t> expected(count);
    std::iota(expected.begin(), expected.end(), 0);
    assert(flattened == expected);
}
}

int main()
{
    const auto blocks = TAP_PBCD::BuildConstantDistanceBlocks(10, 3);
    assert(blocks.size() == 4);
    assert((blocks[0] == std::vector<std::size_t>{0, 3, 6}));
    assert((blocks[1] == std::vector<std::size_t>{1, 4, 7}));
    assert((blocks[2] == std::vector<std::size_t>{2, 5, 8}));
    assert((blocks[3] == std::vector<std::size_t>{9}));
    CheckCoverage(10, 3);
    CheckCoverage(5, 1);
    CheckCoverage(5, 128);
    CheckCoverage(470805, 128);
    assert(TAP_PBCD::BuildConstantDistanceBlocks(0, 128).empty());
    return 0;
}
