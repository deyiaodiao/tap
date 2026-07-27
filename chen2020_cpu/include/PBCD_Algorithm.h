#ifndef CHEN2020_PBCD_ALGORITHM_H
#define CHEN2020_PBCD_ALGORITHM_H

#include "stdafx.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

struct PBCDComponentTiming
{
    double wallSeconds = 0.0;
    double idealEffectiveSeconds = 0.0;
    double serialMergeSeconds = 0.0;
};

struct PBCDTiming
{
    PBCDComponentTiming initialization;
    PBCDComponentTiming pathGeneration;
    PBCDComponentTiming fullAdjustment;
    PBCDComponentTiming restrictedAdjustment;
    PBCDComponentTiming relativeGap;
};

enum class PBCDODScreeningMode
{
    RelativeGap,
    MaximumGPFlowShift
};

class TAP_PBCD : public TNM_TAP
{
public:
    TAP_PBCD();
    ~TAP_PBCD() = default;

    void SetThreadCount(int value);
    void SetGPStep(double value);
    void SetODPerBlock(std::size_t value);
    void SetMaxInnerIterations(int value);
    void SetFullCheckFrequency(int value);
    void SetODScreeningMode(PBCDODScreeningMode value);
    void SetTraceInner(bool value);
    void SetUniformBPR(double alpha, double beta);
    void DisableUniformBPR();

    void PreProcess() override;
    void Initialize() override;
    void MainLoop() override;
    void PostProcess() override;

    double GetRelativeGap() const { return convIndicator; }
    double GetTotalDemand() const { return m_totalDemand; }
    double GetAverageExcessCost() const { return m_averageExcessCost; }
    double GetTotalSystemTravelTime() const { return m_totalSystemTravelTime; }
    double GetShortestPathTotal() const { return m_shortestPathTotal; }
    int GetThreadCount() const { return m_threadCount; }
    std::size_t GetODPerBlock() const { return m_odPerBlock; }
    int GetMaxInnerIterations() const { return m_maxInnerIterations; }
    int GetFullCheckFrequency() const { return m_fullCheckFrequency; }
    double GetGPStep() const { return m_gpStep; }
    PBCDODScreeningMode GetODScreeningMode() const { return m_odScreeningMode; }
    const PBCDTiming& GetTiming() const { return m_timing; }

    bool WriteSolutionJson(const std::string& path, double wallSeconds) const;

    static std::vector<std::vector<std::size_t>> BuildConstantDistanceBlocks(
        std::size_t itemCount,
        std::size_t itemsPerBlock);
    static double ComputeODRelativeGap(
        double demand,
        double shortestPathCost,
        double usedPathTotalCost);

private:
    struct ODRef
    {
        TNM_SORIGIN* origin;
        TNM_SDEST* destination;
    };

    struct ODAdjustment
    {
        std::vector<std::pair<int, double>> linkDeltas;
        double absolutePathFlowChange = 0.0;
        double maximumPathFlowShift = 0.0;
        double odRelativeGap = 0.0;
        bool adjusted = false;
        double taskSeconds = 0.0;
    };

    struct ShortestTree
    {
        std::vector<double> distances;
        std::vector<TNM_SLINK*> predecessors;
    };

    int m_threadCount;
    double m_gpStep;
    std::size_t m_odPerBlock;
    int m_maxInnerIterations;
    int m_fullCheckFrequency;
    PBCDODScreeningMode m_odScreeningMode;
    bool m_traceInner;
    bool m_overrideBPR;
    double m_bprAlpha;
    double m_bprBeta;
    double m_totalDemand;
    double m_totalSystemTravelTime;
    double m_shortestPathTotal;
    double m_averageExcessCost;
    PBCDTiming m_timing;
    std::vector<ODRef> m_odPairs;
    std::vector<std::vector<std::size_t>> m_originODIndices;

    void FilterIntrazonalDemand();
    void ApplyUniformBPR();
    void BuildODIndex();
    bool GeneratePaths(bool initialize, PBCDComponentTiming& timing);
    ShortestTree ComputeShortestPathTree(TNM_SNODE* origin) const;
    std::vector<TNM_SLINK*> ReconstructReversePath(
        TNM_SNODE* origin,
        TNM_SNODE* destination,
        const std::vector<TNM_SLINK*>& predecessors) const;
    bool RouteExists(const TNM_SDEST* destination, const std::vector<TNM_SLINK*>& route) const;
    double PathCost(const TNM_SPATH* path) const;
    double SymmetricDifferenceDerivative(
        const TNM_SPATH* path,
        const TNM_SPATH* shortest) const;
    void AdjustOneOD(
        std::size_t odIndex,
        double maximumODGap,
        ODAdjustment& result);
    double AdjustODSet(
        const std::vector<std::size_t>& odIndices,
        double maximumODGap,
        PBCDComponentTiming& timing,
        std::vector<std::size_t>& restrictedODs);
    void PruneZeroFlowPaths(std::size_t odIndex);
    double ComputeRelativeGapParallel(PBCDComponentTiming& timing);
};

#endif
