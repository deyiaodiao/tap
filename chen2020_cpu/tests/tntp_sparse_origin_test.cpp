#include "PBCD_Algorithm.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace
{
void Require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void WriteFile(const std::filesystem::path& path, const std::string& contents)
{
    std::ofstream stream(path);
    Require(stream.is_open(), "cannot create TNTP test fixture");
    stream << contents;
    Require(stream.good(), "cannot write TNTP test fixture");
}
}

int main()
{
    namespace fs = std::filesystem;
    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    const fs::path directory = fs::temp_directory_path()
        / ("chen2020_sparse_origin_" + std::to_string(nonce));
    fs::create_directories(directory);

    const fs::path prefix = directory / "SparseOrigins";
    WriteFile(prefix.string() + "_net.tntp",
        "<NUMBER OF ZONES> 3\n"
        "<NUMBER OF NODES> 3\n"
        "<FIRST THRU NODE> 1\n"
        "<NUMBER OF LINKS> 3\n"
        "<END OF METADATA>\n"
        "1 2 100 1 1 0.15 4 60 0 1 ;\n"
        "2 3 100 1 1 0.15 4 60 0 1 ;\n"
        "3 1 100 1 1 0.15 4 60 0 1 ;\n");
    WriteFile(prefix.string() + "_trips.tntp",
        "<NUMBER OF ZONES> 3\n"
        "<TOTAL OD FLOW> 12\n"
        "\n"
        "Origin 1\n"
        "3 : 5;\n"
        "Origin 3\n"
        "1 : 7;\n");

    TAP_PBCD solver;
    solver.SetLPF(BPRLK);
    solver.SetCostScalar(1.0);
    solver.SetCostCoef(1.0, 0.0);
    Require(
        solver.Build(prefix.string(), (directory / "result").string(), NETTAPAS) == 0,
        "failed to build sparse-origin TNTP fixture");

    // The metadata declares three zones, but only two real Origin blocks exist.
    // EOF must not create an empty duplicate of the final origin.
    Require(solver.network->numOfOrigin == 2, "EOF created a duplicate final origin");
    Require(solver.network->originVector[0]->id_() == 1, "first origin ID changed");
    Require(solver.network->originVector[1]->id_() == 3, "final origin ID changed");
    Require(solver.network->originVector[1]->numOfDest == 1, "final origin demand was lost");
    Require(
        solver.network->originVector[1]->destVector[0]->dest->id == 1,
        "final origin destination changed");
    Require(
        solver.network->originVector[1]->destVector[0]->assDemand == 7.0,
        "final origin demand changed");

    solver.SetCentroidsBlocked(false);
    solver.PreProcess();
    Require(solver.GetTotalDemand() == 12.0, "preprocessing lost sparse-origin demand");

    fs::remove_all(directory);
    return 0;
}
