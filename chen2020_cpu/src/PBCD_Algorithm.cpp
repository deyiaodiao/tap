#include "stdafx.h"
#include "PBCD_Algorithm.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#else
static double omp_get_wtime()
{
    return static_cast<double>(clock()) / CLOCKS_PER_SEC;
}
static void omp_set_dynamic(int) {}
static void omp_set_num_threads(int) {}
#endif

namespace
{
constexpr double kFlowZeroTolerance = 1e-12;
constexpr double kInnerShiftTolerance = 1e-10;
constexpr std::size_t kDeterministicCommitLanes = 128;

double Maximum(const std::vector<double>& values)
{
    return values.empty() ? 0.0 : *std::max_element(values.begin(), values.end());
}

std::string JsonEscape(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size() + 8);
    for (const char character : value)
    {
        switch (character)
        {
        case '\\': escaped += "\\\\"; break;
        case '"': escaped += "\\\""; break;
        case '\n': escaped += "\\n"; break;
        case '\r': escaped += "\\r"; break;
        case '\t': escaped += "\\t"; break;
        default: escaped += character; break;
        }
    }
    return escaped;
}
}

TAP_PBCD::TAP_PBCD()
    : m_threadCount(1),
      m_gpStep(0.3),
      m_odPerBlock(128),
      m_maxInnerIterations(1000),
      m_fullCheckFrequency(100),
      m_overrideBPR(false),
      m_bprAlpha(0.15),
      m_bprBeta(4.0),
      m_totalDemand(0.0),
      m_totalSystemTravelTime(0.0),
      m_shortestPathTotal(0.0),
      m_averageExcessCost(0.0)
{
}

void TAP_PBCD::SetThreadCount(int value)
{
    if (value < 1) throw std::invalid_argument("thread count must be positive");
    m_threadCount = value;
}

void TAP_PBCD::SetGPStep(double value)
{
    if (!(value > 0.0 && value <= 1.0))
        throw std::invalid_argument("GP step must be in (0, 1]");
    m_gpStep = value;
}

void TAP_PBCD::SetODPerBlock(std::size_t value)
{
    if (value == 0) throw std::invalid_argument("OD per block must be positive");
    m_odPerBlock = value;
}

void TAP_PBCD::SetMaxInnerIterations(int value)
{
    if (value < 0) throw std::invalid_argument("max inner iterations cannot be negative");
    m_maxInnerIterations = value;
}

void TAP_PBCD::SetFullCheckFrequency(int value)
{
    if (value < 1) throw std::invalid_argument("full check frequency must be positive");
    m_fullCheckFrequency = value;
}

void TAP_PBCD::SetUniformBPR(double alpha, double beta)
{
    if (alpha < 0.0 || beta <= 0.0)
        throw std::invalid_argument("invalid uniform BPR parameters");
    m_overrideBPR = true;
    m_bprAlpha = alpha;
    m_bprBeta = beta;
}

void TAP_PBCD::DisableUniformBPR()
{
    m_overrideBPR = false;
}

std::vector<std::vector<std::size_t>> TAP_PBCD::BuildConstantDistanceBlocks(
    std::size_t itemCount,
    std::size_t itemsPerBlock)
{
    if (itemsPerBlock == 0)
        throw std::invalid_argument("itemsPerBlock must be positive");
    std::vector<std::vector<std::size_t>> blocks;
    if (itemCount == 0) return blocks;

    const std::size_t fullBlockCount = itemCount / itemsPerBlock;
    if (fullBlockCount == 0)
    {
        blocks.emplace_back();
        for (std::size_t index = 0; index < itemCount; ++index)
            blocks.back().push_back(index);
        return blocks;
    }

    blocks.reserve(fullBlockCount + 1);
    for (std::size_t block = 0; block < fullBlockCount; ++block)
    {
        blocks.emplace_back();
        blocks.back().reserve(itemsPerBlock);
        for (std::size_t offset = 0; offset < itemsPerBlock; ++offset)
            blocks.back().push_back(block + offset * fullBlockCount);
    }
    const std::size_t firstRemaining = fullBlockCount * itemsPerBlock;
    if (firstRemaining < itemCount)
    {
        blocks.emplace_back();
        for (std::size_t index = firstRemaining; index < itemCount; ++index)
            blocks.back().push_back(index);
    }
    return blocks;
}

void TAP_PBCD::FilterIntrazonalDemand()
{
    for (TNM_SORIGIN* origin : network->originVector)
    {
        double interzonalTotal = 0.0;
        for (int destinationIndex = 0; destinationIndex < origin->numOfDest; ++destinationIndex)
        {
            TNM_SDEST* destination = origin->destVector[destinationIndex];
            if (destination->dest == origin->origin) destination->assDemand = 0.0;
            interzonalTotal += destination->assDemand;
        }
        origin->m_tdmd = interzonalTotal;
    }
    network->ClearZeroDemandOD();
}

void TAP_PBCD::ApplyUniformBPR()
{
    if (!m_overrideBPR) return;
    for (TNM_SLINK* link : network->linkVector)
    {
        TNM_BPRLK* bpr = dynamic_cast<TNM_BPRLK*>(link);
        if (bpr != nullptr) bpr->SetParameters(m_bprAlpha, m_bprBeta);
    }
}

void TAP_PBCD::BuildODIndex()
{
    m_odPairs.clear();
    m_originODIndices.assign(network->numOfOrigin, {});
    m_totalDemand = 0.0;
    for (int originIndex = 0; originIndex < network->numOfOrigin; ++originIndex)
    {
        TNM_SORIGIN* origin = network->originVector[originIndex];
        for (int destinationIndex = 0; destinationIndex < origin->numOfDest; ++destinationIndex)
        {
            TNM_SDEST* destination = origin->destVector[destinationIndex];
            if (destination->assDemand <= 0.0 || destination->dest == origin->origin) continue;
            const std::size_t index = m_odPairs.size();
            m_odPairs.push_back({origin, destination});
            m_originODIndices[originIndex].push_back(index);
            m_totalDemand += destination->assDemand;
        }
    }
}

void TAP_PBCD::PreProcess()
{
    TNM_TAP::PreProcess();
    m_timing = PBCDTiming{};
    omp_set_dynamic(0);
    omp_set_num_threads(m_threadCount);
    network->Reset();
    FilterIntrazonalDemand();
    ApplyUniformBPR();
    network->UpdateLinkCost();
    network->UpdateLinkCostDer();
    BuildODIndex();
}

TAP_PBCD::ShortestTree TAP_PBCD::ComputeShortestPathTree(TNM_SNODE* origin) const
{
    const int nodeCount = network->numOfNode;
    const double infinity = std::numeric_limits<double>::infinity();
    ShortestTree tree;
    tree.distances.assign(nodeCount + 1, infinity);
    tree.predecessors.assign(nodeCount + 1, nullptr);
    using QueueEntry = std::pair<double, int>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    tree.distances[origin->id] = 0.0;
    queue.push({0.0, origin->id});

    while (!queue.empty())
    {
        const QueueEntry current = queue.top();
        queue.pop();
        const double distance = current.first;
        const int nodeId = current.second;
        if (distance != tree.distances[nodeId]) continue;
        TNM_SNODE* node = network->nodeVector[nodeId - 1];
        if (!node->m_isThrough && node != origin) continue;
        for (TNM_SLINK* link : node->forwStar)
        {
            const int head = link->head->id;
            const double candidate = distance + link->cost;
            if (candidate < tree.distances[head])
            {
                tree.distances[head] = candidate;
                tree.predecessors[head] = link;
                queue.push({candidate, head});
            }
        }
    }
    return tree;
}

std::vector<TNM_SLINK*> TAP_PBCD::ReconstructReversePath(
    TNM_SNODE* origin,
    TNM_SNODE* destination,
    const std::vector<TNM_SLINK*>& predecessors) const
{
    std::vector<TNM_SLINK*> route;
    TNM_SNODE* node = destination;
    int guard = 0;
    while (node != origin && guard <= network->numOfNode)
    {
        TNM_SLINK* link = predecessors[node->id];
        if (link == nullptr)
        {
            route.clear();
            return route;
        }
        route.push_back(link);
        node = link->tail;
        ++guard;
    }
    if (node != origin) route.clear();
    return route;
}

bool TAP_PBCD::RouteExists(
    const TNM_SDEST* destination,
    const std::vector<TNM_SLINK*>& route) const
{
    for (const TNM_SPATH* path : destination->pathSet)
    {
        if (path->path.size() != route.size()) continue;
        if (std::equal(path->path.begin(), path->path.end(), route.begin())) return true;
    }
    return false;
}

bool TAP_PBCD::GeneratePaths(bool initialize, PBCDComponentTiming& timing)
{
    const double wallStart = omp_get_wtime();
    if (initialize)
        std::cerr << "STAGE initialization_shortest_paths_begin" << std::endl;
    std::vector<std::vector<TNM_SLINK*>> routes(m_odPairs.size());
    std::vector<double> taskTimes(network->numOfOrigin, 0.0);

#pragma omp parallel for schedule(static) num_threads(m_threadCount)
    for (int originIndex = 0; originIndex < network->numOfOrigin; ++originIndex)
    {
        const double taskStart = omp_get_wtime();
        TNM_SORIGIN* origin = network->originVector[originIndex];
        const ShortestTree tree = ComputeShortestPathTree(origin->origin);
        for (const std::size_t odIndex : m_originODIndices[originIndex])
        {
            routes[odIndex] = ReconstructReversePath(
                origin->origin,
                m_odPairs[odIndex].destination->dest,
                tree.predecessors);
        }
        taskTimes[originIndex] = omp_get_wtime() - taskStart;
    }
    if (initialize)
        std::cerr << "STAGE initialization_shortest_paths_complete wall_seconds="
                  << (omp_get_wtime() - wallStart) << std::endl;

    const double validationStart = omp_get_wtime();
    for (std::size_t odIndex = 0; odIndex < routes.size(); ++odIndex)
    {
        if (routes[odIndex].empty())
        {
            std::cerr << "No path for OD " << m_odPairs[odIndex].origin->origin->id
                      << " -> " << m_odPairs[odIndex].destination->dest->id << std::endl;
            return false;
        }
    }
    const double validationSeconds = omp_get_wtime() - validationStart;

    const std::size_t laneCount = std::min(
        kDeterministicCommitLanes,
        std::max<std::size_t>(1, m_odPairs.size()));
    const std::size_t linkCount = static_cast<std::size_t>(network->numOfLink);
    std::vector<double> laneVolumes;
    if (initialize) laneVolumes.assign(laneCount * linkCount, 0.0);
    std::vector<double> commitTaskTimes(laneCount, 0.0);
    if (initialize)
        std::cerr << "STAGE initialization_path_commit_begin lanes="
                  << laneCount << std::endl;

#pragma omp parallel for schedule(static) num_threads(m_threadCount)
    for (int lane = 0; lane < static_cast<int>(laneCount); ++lane)
    {
        const double taskStart = omp_get_wtime();
        double* localVolumes = initialize
            ? laneVolumes.data() + static_cast<std::size_t>(lane) * linkCount
            : nullptr;
        for (std::size_t odIndex = static_cast<std::size_t>(lane);
             odIndex < m_odPairs.size();
             odIndex += laneCount)
        {
            TNM_SDEST* destination = m_odPairs[odIndex].destination;
            if (!initialize && RouteExists(destination, routes[odIndex])) continue;
            TNM_SPATH* path = new TNM_SPATH;
            path->path = routes[odIndex];
            path->flow = initialize ? destination->assDemand : 0.0;
            destination->pathSet.push_back(path);
            if (initialize)
                for (TNM_SLINK* link : path->path)
                    localVolumes[link->orderID - 1] += path->flow;
        }
        commitTaskTimes[lane] = omp_get_wtime() - taskStart;
    }

    std::vector<double> volumeMergeTaskTimes(initialize ? linkCount : 0, 0.0);
    if (initialize)
    {
#pragma omp parallel for schedule(static) num_threads(m_threadCount)
        for (int linkIndex = 0; linkIndex < static_cast<int>(linkCount); ++linkIndex)
        {
            const double taskStart = omp_get_wtime();
            double volume = 0.0;
            for (std::size_t lane = 0; lane < laneCount; ++lane)
                volume += laneVolumes[lane * linkCount + linkIndex];
            network->linkVector[linkIndex]->volume = volume;
            volumeMergeTaskTimes[linkIndex] = omp_get_wtime() - taskStart;
        }
    }

    if (initialize)
        std::cerr << "STAGE initialization_path_commit_complete wall_seconds="
                  << (omp_get_wtime() - wallStart) << std::endl;

    const double wallSeconds = omp_get_wtime() - wallStart;
    timing.wallSeconds += wallSeconds;
    timing.serialMergeSeconds += validationSeconds;
    timing.idealEffectiveSeconds += Maximum(taskTimes)
        + Maximum(commitTaskTimes)
        + Maximum(volumeMergeTaskTimes)
        + validationSeconds;
    return true;
}

void TAP_PBCD::Initialize()
{
    if (!GeneratePaths(true, m_timing.initialization))
    {
        termFlag = ErrorTerm;
        return;
    }
    network->UpdateLinkCost();
    network->UpdateLinkCostDer();
    ComputeOFV();
    convIndicator = ComputeRelativeGapParallel(m_timing.relativeGap);
    std::cout << std::setprecision(15)
              << "PBCD iteration 0: objective=" << OFV
              << " relative_gap=" << convIndicator
              << " paths=";
    std::size_t pathCount = 0;
    for (const ODRef& od : m_odPairs) pathCount += od.destination->pathSet.size();
    std::cout << pathCount << std::endl;
}

double TAP_PBCD::PathCost(const TNM_SPATH* path) const
{
    double cost = 0.0;
    for (const TNM_SLINK* link : path->path) cost += link->cost;
    return cost;
}

double TAP_PBCD::SymmetricDifferenceDerivative(
    const TNM_SPATH* path,
    const TNM_SPATH* shortest) const
{
    thread_local std::vector<int> linkMarks;
    thread_local int stamp = 0;
    if (linkMarks.size() != static_cast<std::size_t>(network->numOfLink))
    {
        linkMarks.assign(network->numOfLink, 0);
        stamp = 0;
    }
    if (stamp == std::numeric_limits<int>::max())
    {
        std::fill(linkMarks.begin(), linkMarks.end(), 0);
        stamp = 0;
    }
    ++stamp;

    for (const TNM_SLINK* link : path->path)
        linkMarks[link->orderID - 1] = stamp;

    double derivative = 0.0;
    for (const TNM_SLINK* link : shortest->path)
    {
        const int linkIndex = link->orderID - 1;
        if (linkMarks[linkIndex] == stamp)
            linkMarks[linkIndex] = -stamp;
        else
            derivative += link->fdCost;
    }
    for (const TNM_SLINK* link : path->path)
        if (linkMarks[link->orderID - 1] == stamp)
            derivative += link->fdCost;
    return derivative;
}

void TAP_PBCD::AdjustOneOD(
    std::size_t odIndex,
    double maximumODGap,
    ODAdjustment& result)
{
    const double taskStart = omp_get_wtime();
    result.linkDeltas.clear();
    result.absolutePathFlowChange = 0.0;
    result.maximumPathFlowShift = 0.0;
    result.adjusted = false;
    result.taskSeconds = 0.0;
    TNM_SDEST* destination = m_odPairs[odIndex].destination;
    if (destination->pathSet.size() <= 1)
    {
        destination->shiftFlow = 0.0;
        result.taskSeconds = omp_get_wtime() - taskStart;
        return;
    }

    std::size_t shortestIndex = 0;
    double shortestCost = std::numeric_limits<double>::infinity();
    double maximumCost = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < destination->pathSet.size(); ++index)
    {
        TNM_SPATH* path = destination->pathSet[index];
        path->cost = PathCost(path);
        if (path->cost < shortestCost)
        {
            shortestCost = path->cost;
            shortestIndex = index;
        }
        maximumCost = std::max(maximumCost, static_cast<double>(path->cost));
    }
    TNM_SPATH* shortest = destination->pathSet[shortestIndex];

    thread_local std::vector<double> flowReductions;
    flowReductions.assign(destination->pathSet.size(), 0.0);
    for (std::size_t index = 0; index < destination->pathSet.size(); ++index)
    {
        TNM_SPATH* path = destination->pathSet[index];
        if (index == shortestIndex) continue;
        const double gap = std::max(0.0, static_cast<double>(path->cost) - shortestCost);
        const double derivative = SymmetricDifferenceDerivative(path, shortest);
        double flowReduction = 0.0;
        if (gap > 0.0)
        {
            if (derivative > std::numeric_limits<double>::epsilon())
                flowReduction = std::min(
                    static_cast<double>(path->flow),
                    m_gpStep * gap / derivative);
            else
                flowReduction = path->flow;
        }
        flowReductions[index] = flowReduction;
        result.maximumPathFlowShift = std::max(
            result.maximumPathFlowShift,
            flowReduction);
    }
    destination->shiftFlow = result.maximumPathFlowShift;
    if (result.maximumPathFlowShift <= maximumODGap)
    {
        result.taskSeconds = omp_get_wtime() - taskStart;
        return;
    }

    result.adjusted = true;
    std::size_t pathLinkCount = 0;
    for (const TNM_SPATH* path : destination->pathSet)
        pathLinkCount += path->path.size();
    result.linkDeltas.reserve(pathLinkCount);

    double nonShortestTotal = 0.0;
    for (std::size_t index = 0; index < destination->pathSet.size(); ++index)
    {
        TNM_SPATH* path = destination->pathSet[index];
        path->preFlow = path->flow;
        if (index == shortestIndex) continue;
        path->flow = std::max(
            0.0,
            static_cast<double>(path->flow) - flowReductions[index]);
        nonShortestTotal += path->flow;
    }
    shortest->flow = std::max(
        0.0,
        static_cast<double>(destination->assDemand) - nonShortestTotal);

    for (std::size_t index = 0; index < destination->pathSet.size(); ++index)
    {
        TNM_SPATH* path = destination->pathSet[index];
        const double delta = path->flow - path->preFlow;
        result.absolutePathFlowChange += std::abs(delta);
        if (std::abs(delta) <= kFlowZeroTolerance) continue;
        for (TNM_SLINK* link : path->path)
            result.linkDeltas.push_back({link->orderID, delta});
    }
    result.taskSeconds = omp_get_wtime() - taskStart;
}

void TAP_PBCD::PruneZeroFlowPaths(std::size_t odIndex)
{
    TNM_SDEST* destination = m_odPairs[odIndex].destination;
    auto iterator = destination->pathSet.begin();
    while (iterator != destination->pathSet.end() && destination->pathSet.size() > 1)
    {
        if ((*iterator)->flow <= kFlowZeroTolerance)
        {
            delete *iterator;
            iterator = destination->pathSet.erase(iterator);
        }
        else ++iterator;
    }
}

double TAP_PBCD::AdjustODSet(
    const std::vector<std::size_t>& odIndices,
    double maximumODGap,
    PBCDComponentTiming& timing,
    std::vector<std::size_t>& restrictedODs)
{
    restrictedODs.clear();
    if (restrictedODs.capacity() < odIndices.size())
        restrictedODs.reserve(odIndices.size());
    if (odIndices.empty()) return 0.0;
    const double wallStart = omp_get_wtime();
    double idealSeconds = 0.0;
    double serialSeconds = 0.0;
    double totalShift = 0.0;
    const std::size_t itemCount = odIndices.size();
    const std::size_t fullBlockCount = itemCount / m_odPerBlock;
    const std::size_t remainingCount = fullBlockCount == 0
        ? itemCount
        : itemCount - fullBlockCount * m_odPerBlock;
    const std::size_t blockCount = fullBlockCount == 0
        ? 1
        : fullBlockCount + (remainingCount > 0 ? 1 : 0);
    std::vector<int> touchedStamp(network->numOfLink, -1);
    int blockStamp = 0;

    const std::size_t maximumBlockSize = std::min(itemCount, m_odPerBlock);
    std::vector<ODAdjustment> adjustments(maximumBlockSize);
    std::vector<std::size_t> committedODs(maximumBlockSize);
    std::vector<double> postProcessTimes(maximumBlockSize, 0.0);
    std::vector<int> touchedLinks;
    touchedLinks.reserve(network->numOfLink);
    std::vector<double> linkUpdateTimes(network->numOfLink, 0.0);
    double longestTask = 0.0;
    double mergeSeconds = 0.0;
    bool abortBlocks = false;

#pragma omp parallel num_threads(m_threadCount) shared(abortBlocks, blockStamp)
    {
        for (std::size_t blockIndex = 0; blockIndex < blockCount; ++blockIndex)
        {
            const bool remainingBlock = fullBlockCount > 0
                && blockIndex == fullBlockCount;
            const std::size_t blockSize = fullBlockCount == 0
                ? itemCount
                : (remainingBlock ? remainingCount : m_odPerBlock);
#pragma omp for schedule(static)
            for (int position = 0; position < static_cast<int>(blockSize); ++position)
            {
                const std::size_t inputPosition = fullBlockCount == 0 || remainingBlock
                    ? (remainingBlock
                        ? fullBlockCount * m_odPerBlock + position
                        : position)
                    : blockIndex + static_cast<std::size_t>(position) * fullBlockCount;
                AdjustOneOD(
                    odIndices[inputPosition],
                    maximumODGap,
                    adjustments[position]);
            }

#pragma omp single
            {
                longestTask = 0.0;
                for (std::size_t position = 0; position < blockSize; ++position)
                    longestTask = std::max(longestTask, adjustments[position].taskSeconds);

                touchedLinks.clear();
                const double mergeStart = omp_get_wtime();
                for (std::size_t position = 0; position < blockSize; ++position)
                {
                    const std::size_t inputPosition = fullBlockCount == 0 || remainingBlock
                        ? (remainingBlock
                            ? fullBlockCount * m_odPerBlock + position
                            : position)
                        : blockIndex + position * fullBlockCount;
                    const std::size_t odIndex = odIndices[inputPosition];
                    committedODs[position] = odIndex;
                    const ODAdjustment& adjustment = adjustments[position];
                    totalShift += adjustment.absolutePathFlowChange;
                    if (adjustment.adjusted)
                        restrictedODs.push_back(odIndex);
                    for (const auto& delta : adjustment.linkDeltas)
                    {
                        TNM_SLINK* link = network->linkVector[delta.first - 1];
                        link->volume += delta.second;
                        if (std::abs(link->volume) <= kFlowZeroTolerance) link->volume = 0.0;
                        if (link->volume < -1e-8)
                        {
                            std::cerr << "Negative link volume after PBCD merge on link "
                                      << link->id << ": " << link->volume << std::endl;
                            termFlag = ErrorTerm;
                        }
                        const int linkIndex = delta.first - 1;
                        if (touchedStamp[linkIndex] != blockStamp)
                        {
                            touchedStamp[linkIndex] = blockStamp;
                            touchedLinks.push_back(linkIndex);
                        }
                    }
                }
                mergeSeconds = omp_get_wtime() - mergeStart;
            }

#pragma omp for schedule(static)
            for (int position = 0; position < static_cast<int>(touchedLinks.size()); ++position)
            {
                const double taskStart = omp_get_wtime();
                const int linkIndex = touchedLinks[position];
                TNM_SLINK* link = network->linkVector[linkIndex];
                link->cost = link->GetCost();
                link->fdCost = link->GetDerCost();
                linkUpdateTimes[position] = omp_get_wtime() - taskStart;
            }

#pragma omp for schedule(static)
            for (int position = 0; position < static_cast<int>(blockSize); ++position)
            {
                const double taskStart = omp_get_wtime();
                const std::size_t odIndex = committedODs[position];
                if (adjustments[position].adjusted)
                    PruneZeroFlowPaths(odIndex);
                postProcessTimes[position] = omp_get_wtime() - taskStart;
            }

#pragma omp single
            {
                double longestLinkUpdate = 0.0;
                for (std::size_t position = 0; position < touchedLinks.size(); ++position)
                    longestLinkUpdate = std::max(longestLinkUpdate, linkUpdateTimes[position]);
                double longestPostProcess = 0.0;
                for (std::size_t position = 0; position < blockSize; ++position)
                    longestPostProcess = std::max(
                        longestPostProcess, postProcessTimes[position]);
                serialSeconds += mergeSeconds;
                idealSeconds += longestTask + mergeSeconds
                    + longestLinkUpdate + longestPostProcess;
                abortBlocks = termFlag == ErrorTerm;
                if (!abortBlocks) ++blockStamp;
            }
            if (abortBlocks) break;
        }
    }

    timing.wallSeconds += omp_get_wtime() - wallStart;
    timing.idealEffectiveSeconds += idealSeconds;
    timing.serialMergeSeconds += serialSeconds;
    return totalShift;
}

double TAP_PBCD::ComputeRelativeGapParallel(PBCDComponentTiming& timing)
{
    const double wallStart = omp_get_wtime();
    std::vector<double> originTerms(network->numOfOrigin, 0.0);
    std::vector<double> taskTimes(network->numOfOrigin, 0.0);

#pragma omp parallel for schedule(static) num_threads(m_threadCount)
    for (int originIndex = 0; originIndex < network->numOfOrigin; ++originIndex)
    {
        const double taskStart = omp_get_wtime();
        TNM_SORIGIN* origin = network->originVector[originIndex];
        const ShortestTree tree = ComputeShortestPathTree(origin->origin);
        double term = 0.0;
        for (const std::size_t odIndex : m_originODIndices[originIndex])
        {
            TNM_SDEST* destination = m_odPairs[odIndex].destination;
            const double pathCost = tree.distances[destination->dest->id];
            if (!std::isfinite(pathCost))
            {
                term = std::numeric_limits<double>::infinity();
                break;
            }
            term += destination->assDemand * pathCost;
        }
        originTerms[originIndex] = term;
        taskTimes[originIndex] = omp_get_wtime() - taskStart;
    }

    const double serialStart = omp_get_wtime();
    m_shortestPathTotal = std::accumulate(originTerms.begin(), originTerms.end(), 0.0);
    m_totalSystemTravelTime = 0.0;
    for (TNM_SLINK* link : network->linkVector)
        m_totalSystemTravelTime += link->volume * link->cost;
    const double excess = m_totalSystemTravelTime - m_shortestPathTotal;
    m_averageExcessCost = m_totalDemand > 0.0 ? excess / m_totalDemand : 0.0;
    const double signedGap = m_totalSystemTravelTime > 0.0
        ? excess / m_totalSystemTravelTime
        : 0.0;
    const double serialSeconds = omp_get_wtime() - serialStart;
    timing.wallSeconds += omp_get_wtime() - wallStart;
    timing.serialMergeSeconds += serialSeconds;
    timing.idealEffectiveSeconds += Maximum(taskTimes) + serialSeconds;
    return std::abs(signedGap);
}

void TAP_PBCD::MainLoop()
{
    const double previousGap = convIndicator;
    const double pathGenerationBefore = m_timing.pathGeneration.wallSeconds;
    const double fullAdjustmentBefore = m_timing.fullAdjustment.wallSeconds;
    const double restrictedAdjustmentBefore = m_timing.restrictedAdjustment.wallSeconds;
    const double relativeGapBefore = m_timing.relativeGap.wallSeconds;
    if (!GeneratePaths(false, m_timing.pathGeneration))
    {
        termFlag = ErrorTerm;
        return;
    }

    std::vector<std::size_t> allODs(m_odPairs.size());
    std::iota(allODs.begin(), allODs.end(), 0);
    const double selectionThreshold = previousGap / 2.0;
    int innerIterationsExecuted = 0;
    std::size_t maximumRestrictedODs = 0;
    std::vector<std::size_t> restricted;
    std::vector<std::size_t> nextRestricted;
    for (int inner = 0; inner < m_maxInnerIterations; ++inner)
    {
        const bool fullCheck = inner % m_fullCheckFrequency == 0;
        const std::vector<std::size_t>& inputODs = fullCheck ? allODs : restricted;
        if (inputODs.empty()) break;
        const double shifted = AdjustODSet(
            inputODs,
            selectionThreshold,
            fullCheck ? m_timing.fullAdjustment : m_timing.restrictedAdjustment,
            nextRestricted);
        if (termFlag == ErrorTerm) return;
        ++innerIterationsExecuted;
        maximumRestrictedODs = std::max(
            maximumRestrictedODs,
            nextRestricted.size());
        restricted.swap(nextRestricted);
        if (restricted.empty()
            || shifted < kInnerShiftTolerance)
            break;
    }

    convIndicator = ComputeRelativeGapParallel(m_timing.relativeGap);
    ComputeOFV();
    std::size_t pathCount = 0;
    for (const ODRef& od : m_odPairs) pathCount += od.destination->pathSet.size();
    std::cout << std::setprecision(15)
              << "PBCD iteration " << curIter
              << ": objective=" << OFV
              << " relative_gap=" << convIndicator
              << " aec=" << m_averageExcessCost
              << " paths=" << pathCount << std::endl;
    std::cerr << std::setprecision(10)
              << "PBCD_COMPONENT iteration=" << curIter
              << " path_generation="
              << (m_timing.pathGeneration.wallSeconds - pathGenerationBefore)
              << " full_adjustment="
              << (m_timing.fullAdjustment.wallSeconds - fullAdjustmentBefore)
              << " restricted_adjustment="
              << (m_timing.restrictedAdjustment.wallSeconds - restrictedAdjustmentBefore)
              << " relative_gap="
              << (m_timing.relativeGap.wallSeconds - relativeGapBefore)
              << " inner_iterations=" << innerIterationsExecuted
              << " max_restricted_ods=" << maximumRestrictedODs << std::endl;
}

void TAP_PBCD::PostProcess()
{
    network->UpdateLinkCost();
    network->UpdateLinkCostDer();
}

bool TAP_PBCD::WriteSolutionJson(const std::string& path, double wallSeconds) const
{
    std::ofstream output(path.c_str());
    if (!output.is_open()) return false;
    output << std::setprecision(17);
    output << "{\n"
           << "  \"instance\": \"" << JsonEscape(inFileName) << "\",\n"
           << "  \"method\": \"chen2020_pbcd\",\n"
           << "  \"implementation\": {\"version\": \"chen2020_cpu-pbcd-v2\", "
           << "\"compiler\": \"" << JsonEscape(__VERSION__) << "\", "
           << "\"cplusplus\": " << __cplusplus
#ifdef _OPENMP
           << ", \"openmp\": " << _OPENMP
#else
           << ", \"openmp\": null"
#endif
           << "},\n"
           << "  \"status\": \""
           << (termFlag == ConvergeTerm ? "completed" :
               termFlag == MaxIterTerm ? "max_iterations" : "error") << "\",\n"
           << "  \"objective_reported\": " << OFV << ",\n"
           << "  \"relative_gap_reported\": " << convIndicator << ",\n"
           << "  \"od_flows\": [\n";

    for (std::size_t odIndex = 0; odIndex < m_odPairs.size(); ++odIndex)
    {
        const ODRef& od = m_odPairs[odIndex];
        output << "    {\"origin\": " << od.origin->origin->id
               << ", \"destination\": " << od.destination->dest->id
               << ", \"demand\": " << od.destination->assDemand
               << ", \"paths\": [";
        for (std::size_t pathIndex = 0; pathIndex < od.destination->pathSet.size(); ++pathIndex)
        {
            const TNM_SPATH* route = od.destination->pathSet[pathIndex];
            if (pathIndex != 0) output << ", ";
            output << "{\"flow\": " << route->flow << ", \"links\": [";
            for (auto iterator = route->path.rbegin(); iterator != route->path.rend(); ++iterator)
            {
                if (iterator != route->path.rbegin()) output << ", ";
                // Serialize the stable TNTP data-row position, not TNM's internal
                // id. TNM also counts <ORIGINAL HEADER> as an id slot on some
                // networks (notably Winnipeg), so internal ids are not portable.
                output << (*iterator)->orderID;
            }
            output << "]}";
        }
        output << "]}" << (odIndex + 1 == m_odPairs.size() ? "\n" : ",\n");
    }

    const auto writeTiming = [&output](const char* name, const PBCDComponentTiming& value, bool comma)
    {
        output << "    \"" << name << "\": {\"wall_seconds\": "
               << value.wallSeconds << ", \"ideal_effective_seconds\": "
               << value.idealEffectiveSeconds << ", \"serial_merge_seconds\": "
               << value.serialMergeSeconds << "}" << (comma ? "," : "") << "\n";
    };

    output << "  ],\n"
           << "  \"configuration\": {\"threads\": " << m_threadCount
           << ", \"convergence_target\": " << convCriterion
           << ", \"max_outer_iterations\": " << maxMainIter
           << ", \"outer_iterations_completed\": " << curIter
           << ", \"cost_scalar\": " << costScalar
           << ", \"seed\": null"
           << ", \"gp_step\": " << m_gpStep
           << ", \"od_per_block\": " << m_odPerBlock
           << ", \"max_inner_iterations\": " << m_maxInnerIterations
           << ", \"full_check_frequency\": " << m_fullCheckFrequency
           << ", \"restricted_od_metric\": \"maximum_gp_path_flow_shift\""
           << ", \"uniform_bpr\": " << (m_overrideBPR ? "true" : "false")
           << ", \"bpr_alpha\": " << m_bprAlpha
           << ", \"bpr_beta\": " << m_bprBeta << "},\n"
           << "  \"metrics\": {\"total_demand\": " << m_totalDemand
           << ", \"total_system_travel_time\": " << m_totalSystemTravelTime
           << ", \"shortest_path_total\": " << m_shortestPathTotal
           << ", \"average_excess_cost\": " << m_averageExcessCost << "},\n"
           << "  \"timing\": {\n";
    writeTiming("initialization", m_timing.initialization, true);
    writeTiming("path_generation", m_timing.pathGeneration, true);
    writeTiming("full_adjustment", m_timing.fullAdjustment, true);
    writeTiming("restricted_adjustment", m_timing.restrictedAdjustment, true);
    writeTiming("relative_gap", m_timing.relativeGap, true);
    output << "    \"wall_seconds\": " << wallSeconds << "\n"
           << "  }\n"
           << "}\n";
    return output.good();
}
