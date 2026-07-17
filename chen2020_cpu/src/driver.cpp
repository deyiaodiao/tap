#include "stdafx.h"
#include "PBCD_Algorithm.h"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#else
static double omp_get_wtime()
{
    return static_cast<double>(clock()) / CLOCKS_PER_SEC;
}
#endif

namespace fs = std::filesystem;

namespace
{
struct Options
{
    int threads = 1;
    double convergence = 1e-12;
    int maxOuter = 500;
    double gpStep = 0.3;
    std::size_t odPerBlock = 128;
    int maxInner = 1000;
    int fullCheckFrequency = 100;
    double costScalar = 60.0;
    bool alphaSet = false;
    bool betaSet = false;
    double bprAlpha = 0.15;
    double bprBeta = 4.0;
    bool loadOnly = false;
    bool initializeOnly = false;
    bool skipArtifacts = false;
    bool skipSolutionJson = false;
};

const char* TerminationName(TERMFLAGS term)
{
    switch (term)
    {
    case ConvergeTerm: return "converged";
    case MaxIterTerm: return "max_iterations";
    case UserTerm: return "user";
    case ErrorTerm: return "error";
    default: return "other";
    }
}

void Usage(const char* executable)
{
    std::cout
        << "Usage: " << executable
        << " pbcd <network-directory> <network-name> <output-prefix> [options]\n"
        << "Options:\n"
        << "  --threads N\n"
        << "  --convergence X\n"
        << "  --max-outer N\n"
        << "  --gp-step X\n"
        << "  --od-per-block N\n"
        << "  --max-inner N\n"
        << "  --full-check-frequency N\n"
        << "  --cost-scalar X\n"
        << "  --load-only 0|1\n"
        << "  --initialize-only 0|1\n"
        << "  --skip-artifacts 0|1\n"
        << "  --skip-solution-json 0|1\n"
        << "  --bpr-alpha X --bpr-beta X\n";
}

bool IsWithin(const fs::path& candidate, const fs::path& root)
{
    const fs::path relative = candidate.lexically_relative(root);
    return !relative.empty() && *relative.begin() != "..";
}

Options ParseOptions(int argc, char* argv[], int first)
{
    Options options;
    for (int index = first; index < argc; index += 2)
    {
        if (index + 1 >= argc)
            throw std::invalid_argument(std::string("missing value for ") + argv[index]);
        const std::string name = argv[index];
        const std::string value = argv[index + 1];
        if (name == "--threads") options.threads = std::stoi(value);
        else if (name == "--convergence") options.convergence = std::stod(value);
        else if (name == "--max-outer") options.maxOuter = std::stoi(value);
        else if (name == "--gp-step") options.gpStep = std::stod(value);
        else if (name == "--od-per-block") options.odPerBlock = std::stoull(value);
        else if (name == "--max-inner") options.maxInner = std::stoi(value);
        else if (name == "--full-check-frequency")
            options.fullCheckFrequency = std::stoi(value);
        else if (name == "--cost-scalar") options.costScalar = std::stod(value);
        else if (name == "--load-only") options.loadOnly = std::stoi(value) != 0;
        else if (name == "--initialize-only")
            options.initializeOnly = std::stoi(value) != 0;
        else if (name == "--skip-artifacts")
            options.skipArtifacts = std::stoi(value) != 0;
        else if (name == "--skip-solution-json")
            options.skipSolutionJson = std::stoi(value) != 0;
        else if (name == "--bpr-alpha")
        {
            options.bprAlpha = std::stod(value);
            options.alphaSet = true;
        }
        else if (name == "--bpr-beta")
        {
            options.bprBeta = std::stod(value);
            options.betaSet = true;
        }
        else throw std::invalid_argument("unknown option: " + name);
    }
    if (options.alphaSet != options.betaSet)
        throw std::invalid_argument("--bpr-alpha and --bpr-beta must be supplied together");
    return options;
}
}

int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--help")
    {
        Usage(argv[0]);
        return 0;
    }
    if (argc < 5 || std::string(argv[1]) != "pbcd")
    {
        Usage(argv[0]);
        return 1;
    }

    try
    {
        std::string networkDirectory = argv[2];
        if (!networkDirectory.empty() && networkDirectory.back() != '/')
            networkDirectory += '/';
        const std::string networkName = argv[3];
        const Options options = ParseOptions(argc, argv, 5);

        const fs::path inputDirectory = fs::weakly_canonical(fs::path(networkDirectory));
        fs::path frozenRoot = inputDirectory;
        if (inputDirectory.parent_path().filename() == "Network")
            frozenRoot = inputDirectory.parent_path().parent_path();
        const fs::path outputPrefix = fs::absolute(fs::path(argv[4])).lexically_normal();
        if (IsWithin(outputPrefix, frozenRoot))
        {
            std::cerr << "Refusing to write output under frozen input tree: "
                      << frozenRoot << std::endl;
            return 4;
        }
        if (!outputPrefix.parent_path().empty())
            fs::create_directories(outputPrefix.parent_path());

        TNM_FloatFormat::SetFormat(18, 10);
        TAP_PBCD solver;
        solver.SetConv(options.convergence);
        solver.SetMaxIter(options.maxOuter);
        solver.SetLPF(BPRLK);
        solver.SetCostScalar(options.costScalar);
        solver.SetCostCoef(1.0, 0.0);
        solver.SetThreadCount(options.threads);
        solver.SetGPStep(options.gpStep);
        solver.SetODPerBlock(options.odPerBlock);
        solver.SetMaxInnerIterations(options.maxInner);
        solver.SetFullCheckFrequency(options.fullCheckFrequency);
        if (options.alphaSet)
            solver.SetUniformBPR(options.bprAlpha, options.bprBeta);

        const std::string inputPrefix = networkDirectory + networkName;
        std::cerr << "STAGE build_begin" << std::endl;
        const double buildStart = omp_get_wtime();
        if (solver.Build(inputPrefix, outputPrefix.string(), NETTAPAS) != 0)
        {
            std::cerr << "Failed to build network object" << std::endl;
            return 2;
        }
        std::cerr << "STAGE build_complete wall_seconds="
                  << std::setprecision(17) << (omp_get_wtime() - buildStart)
                  << std::endl;
        if (options.loadOnly)
        {
            std::cerr << "STAGE load_only_complete nodes=" << solver.network->numOfNode
                      << " links=" << solver.network->numOfLink
                      << " origins=" << solver.network->numOfOrigin << std::endl;
            return 0;
        }
        solver.SetTollType(TT_NOTOLL);
        solver.SetCentroidsBlocked(false);

        if (options.initializeOnly)
        {
            const double initializationStart = omp_get_wtime();
            std::cerr << "STAGE initialize_only_begin" << std::endl;
            solver.PreProcess();
            solver.Initialize();
            solver.PostProcess();
            std::cerr << "STAGE initialize_only_complete wall_seconds="
                      << std::setprecision(17)
                      << (omp_get_wtime() - initializationStart) << std::endl;
            return solver.termFlag == ErrorTerm ? 3 : 0;
        }

        const double wallStart = omp_get_wtime();
        std::cerr << "STAGE solve_begin" << std::endl;
        const TERMFLAGS term = solver.Solve();
        const double wallSeconds = omp_get_wtime() - wallStart;
        std::cerr << "STAGE solve_complete wall_seconds="
                  << std::setprecision(17) << wallSeconds << std::endl;
        if (options.skipArtifacts)
        {
            const PBCDTiming& timing = solver.GetTiming();
            const auto printTiming = [](const char* name, const PBCDComponentTiming& value)
            {
                std::cout << "TIMING_" << name << "_WALL=" << value.wallSeconds << '\n'
                          << "TIMING_" << name << "_IDEAL=" << value.idealEffectiveSeconds << '\n'
                          << "TIMING_" << name << "_SERIAL=" << value.serialMergeSeconds << '\n';
            };
            std::cout << std::setprecision(17)
                      << "OFV=" << solver.OFV << '\n'
                      << "RELATIVE_GAP=" << solver.GetRelativeGap() << '\n'
                      << "ITERATIONS=" << solver.curIter << '\n'
                      << "WALL_SECONDS=" << wallSeconds << '\n';
            printTiming("INITIALIZATION", timing.initialization);
            printTiming("PATH_GENERATION", timing.pathGeneration);
            printTiming("FULL_ADJUSTMENT", timing.fullAdjustment);
            printTiming("RESTRICTED_ADJUSTMENT", timing.restrictedAdjustment);
            printTiming("RELATIVE_GAP", timing.relativeGap);
            std::cout << "TERMINATION=" << TerminationName(term) << std::endl;
            return term == ErrorTerm ? 3 : 0;
        }
        if (!options.skipSolutionJson
            && !solver.WriteSolutionJson(outputPrefix.string() + ".json", wallSeconds))
        {
            std::cerr << "Failed to write solution JSON" << std::endl;
            return 6;
        }
        solver.reportIterHistory = true;
        solver.reportLinkDetail = true;
        solver.reportPathDetail = false;
        solver.Report();

        std::cout << std::setprecision(17)
                  << "OFV=" << solver.OFV << '\n'
                  << "RELATIVE_GAP=" << solver.GetRelativeGap() << '\n'
                  << "AVERAGE_EXCESS_COST=" << solver.GetAverageExcessCost() << '\n'
                  << "TOTAL_DEMAND=" << solver.GetTotalDemand() << '\n'
                  << "ITERATIONS=" << solver.curIter << '\n'
                  << "CPU_SECONDS=" << solver.cpuTime << '\n'
                  << "WALL_SECONDS=" << wallSeconds << '\n'
                  << "THREADS=" << solver.GetThreadCount() << '\n'
                  << "OD_PER_BLOCK=" << solver.GetODPerBlock() << '\n'
                  << "TERMINATION=" << TerminationName(term) << std::endl;
        if (term == ErrorTerm) return 3;
        return term == ConvergeTerm ? 0 : 5;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Invalid configuration or runtime error: " << error.what() << std::endl;
        return 1;
    }
}
