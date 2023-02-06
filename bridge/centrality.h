
#ifndef NK_CENTRALITY_H
#define NK_CENTRALITY_H

#include "rust/cxx.h"
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/centrality/ApproxGroupBetweenness.hpp>
#include <networkit/centrality/ApproxSpanningEdge.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/centrality/DegreeCentrality.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/DynBetweenness.hpp>
#include <networkit/centrality/DynBetweennessOneNode.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>
#include <networkit/centrality/DynTopHarmonicCloseness.hpp>
#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/ForestCentrality.hpp>
#include <networkit/centrality/GedWalk.hpp>
#include "graph_event.h"

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<ApproxBetweenness> NewApproxBetweenness(
        const Graph &G, double epsilon = 0.01, double delta = 0.1,
        double universalConstant = 1.0)
    {
        return make_unique<ApproxBetweenness>(G, epsilon, delta, universalConstant);
    }
    inline void ApproxBetweennessRanking(ApproxBetweenness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ApproxBetweennessScores(ApproxBetweenness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<ApproxCloseness> NewApproxCloseness(
        const Graph &G, count nSamples, double epsilon, bool normalized,
        uint8_t t)
    {
        ApproxCloseness::ClosenessType type;
        switch (t)
        {
        case 0:
            type = ApproxCloseness::ClosenessType::OUTBOUND;
            break;
        case 1:
            type = ApproxCloseness::ClosenessType::INBOUND;
            break;
        case 2:
            type = ApproxCloseness::ClosenessType::SUM;
            break;
        default:
            type = ApproxCloseness::ClosenessType::OUTBOUND;
        }
        return make_unique<ApproxCloseness>(G, nSamples, epsilon, normalized, type);
    }
    inline void ApproxClosenessRanking(ApproxCloseness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ApproxClosenessScores(ApproxCloseness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<double>> ApproxClosenessGetSquareErrorEstimates(ApproxCloseness &algo)
    {
        return make_unique<vector<double>>(algo.getSquareErrorEstimates());
    }

    inline unique_ptr<ApproxElectricalCloseness> NewApproxElectricalCloseness(
        const Graph &G, double epsilon = 0.1, double kappa = 0.3)
    {
        return make_unique<ApproxElectricalCloseness>(G, epsilon, kappa);
    }
    inline void ApproxElectricalClosenessRanking(ApproxElectricalCloseness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ApproxElectricalClosenessScores(ApproxElectricalCloseness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<double>> ApproxElectricalClosenessComputeExactDiagonal(const ApproxElectricalCloseness &algo, double tol)
    {
        return make_unique<vector<double>>(algo.computeExactDiagonal(tol));
    }

    inline unique_ptr<vector<double>> ApproxElectricalClosenessGetDiagonal(const ApproxElectricalCloseness &algo)
    {
        return make_unique<vector<double>>(algo.getDiagonal());
    }

    inline unique_ptr<ApproxGroupBetweenness> NewApproxGroupBetweenness(
        const Graph &G, count groupSize, double epsilon)
    {
        return make_unique<ApproxGroupBetweenness>(G, groupSize, epsilon);
    }

    inline unique_ptr<vector<node>> ApproxGroupBetweennessGroupMaxBetweenness(const ApproxGroupBetweenness &algo)
    {
        return make_unique<vector<node>>(algo.groupMaxBetweenness());
    }

    inline unique_ptr<vector<node>> ApproxGroupBetweennessScoreOfGroup(const ApproxGroupBetweenness &algo,
                                                                       rust::Slice<const node> S, bool normalized)
    {
        vector<node> nodes(S.begin(), S.end());
        return make_unique<vector<node>>(algo.scoreOfGroup(nodes, normalized));
    }

    inline unique_ptr<ApproxSpanningEdge> NewApproxSpanningEdge(
        const Graph &G, double epsilon)
    {
        return make_unique<ApproxSpanningEdge>(G, epsilon);
    }

    inline unique_ptr<vector<double>> ApproxSpanningEdgeScores(const ApproxSpanningEdge &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<Betweenness> NewBetweenness(
        const Graph &G, bool normalized = false, bool computeEdgeCentrality = false)
    {
        return make_unique<Betweenness>(G, normalized, computeEdgeCentrality);
    }
    inline void BetweennessRanking(Betweenness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> BetweennessScores(Betweenness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<double>> BetweennessEdgeScores(Betweenness &algo)
    {
        return make_unique<vector<double>>(algo.edgeScores());
    }

    inline unique_ptr<Closeness> NewCloseness(
        const Graph &G, bool normalized, uint8_t variant)
    {
        return make_unique<Closeness>(G, normalized, variant);
    }
    inline void ClosenessRanking(Closeness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ClosenessScores(Closeness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<CoreDecomposition> NewCoreDecomposition(
        const Graph &G, bool normalized = false,
        bool enforceBucketQueueAlgorithm = false, bool storeNodeOrder = false)
    {
        return make_unique<CoreDecomposition>(G, normalized, enforceBucketQueueAlgorithm, storeNodeOrder);
    }
    inline void CoreDecompositionRanking(CoreDecomposition &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> CoreDecompositionScores(CoreDecomposition &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<Cover> CoreDecompositionGetCover(const CoreDecomposition &algo)
    {
        return make_unique<Cover>(algo.getCover());
    }

    inline unique_ptr<vector<node>> CoreDecompositionGetNodeOrder(const CoreDecomposition &algo)
    {
        return make_unique<vector<node>>(algo.getNodeOrder());
    }

    inline unique_ptr<Partition> CoreDecompositionGetPartition(const CoreDecomposition &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline unique_ptr<DegreeCentrality> NewDegreeCentrality(
        const Graph &G, bool normalized = false,
        bool outDeg = false, bool ignoreSelfLoops = false)
    {
        return make_unique<DegreeCentrality>(G, normalized, outDeg, ignoreSelfLoops);
    }
    inline void DegreeCentralityRanking(DegreeCentrality &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> DegreeCentralityScores(DegreeCentrality &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<DynApproxBetweenness> NewDynApproxBetweenness(
        const Graph &G, double epsilon = 0.01, double delta = 0.1,
        bool storePredecessors = true, double universalConstant = 0.5)
    {
        return make_unique<DynApproxBetweenness>(G, epsilon, delta, storePredecessors, universalConstant);
    }

    inline void DynApproxBetweennessRanking(DynApproxBetweenness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> DynApproxBetweennessScores(DynApproxBetweenness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline void DynApproxBetweennessUpdate(DynApproxBetweenness &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynApproxBetweennessUpdateBatch(DynApproxBetweenness &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<DynBetweenness> NewDynBetweenness(
        const Graph &G)
    {
        return make_unique<DynBetweenness>(G);
    }

    inline void DynBetweennessRanking(DynBetweenness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> DynBetweennessScores(DynBetweenness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline void DynBetweennessUpdate(DynBetweenness &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynBetweennessUpdateBatch(DynBetweenness &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<DynBetweennessOneNode> NewDynBetweennessOneNode(
        Graph &G, node x)
    {
        return make_unique<DynBetweennessOneNode>(G, x);
    }

    inline void DynBetweennessOneNodeUpdate(DynBetweennessOneNode &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline edgeweight DynBetweennessOneNodeComputeScore(DynBetweennessOneNode &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        return algo.computeScore(toGraphEvent(kind, u, v, ew));
    }

    inline void DynBetweennessOneNodeUpdateBatch(DynBetweennessOneNode &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<DynKatzCentrality> NewDynKatzCentrality(
        const Graph &G, count k, bool groupOnly = false, double tolerance = 1e-9)
    {
        return make_unique<DynKatzCentrality>(G, k, groupOnly, tolerance);
    }

    inline void DynKatzCentralityRanking(DynKatzCentrality &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> DynKatzCentralityScores(DynKatzCentrality &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<node>> DynKatzCentralityTop(DynKatzCentrality &algo, count n)
    {
        return make_unique<vector<node>>(algo.top(n));
    }

    inline void DynKatzCentralityUpdate(DynKatzCentrality &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynKatzCentralityUpdateBatch(DynKatzCentrality &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<DynTopHarmonicCloseness> NewDynTopHarmonicCloseness(
        const Graph &G, count k = 1, bool useBFSbound = false)
    {
        return make_unique<DynTopHarmonicCloseness>(G, k, useBFSbound);
    }

    inline void DynTopHarmonicClosenessRanking(DynTopHarmonicCloseness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline void DynTopHarmonicClosenessUpdate(DynTopHarmonicCloseness &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynTopHarmonicClosenessUpdateBatch(DynTopHarmonicCloseness &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<vector<node>> DynTopHarmonicClosenessTopkNodesList(DynTopHarmonicCloseness &algo, bool includeTrail)
    {
        return make_unique<vector<node>>(algo.topkNodesList(includeTrail));
    }

    inline unique_ptr<vector<edgeweight>> DynTopHarmonicClosenessTopkScoresList(DynTopHarmonicCloseness &algo, bool includeTrail)
    {
        return make_unique<vector<edgeweight>>(algo.topkScoresList(includeTrail));
    }

    inline unique_ptr<EigenvectorCentrality> NewEigenvectorCentrality(
        const Graph &G, double tol)
    {
        return make_unique<EigenvectorCentrality>(G, tol);
    }
    inline void EigenvectorCentralityRanking(EigenvectorCentrality &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> EigenvectorCentralityScores(EigenvectorCentrality &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<EstimateBetweenness> NewEstimateBetweenness(
        const Graph &G, count nSamples, bool normalized = false,
        bool parallel_flag = false)
    {
        return make_unique<EstimateBetweenness>(G, nSamples, normalized, parallel_flag);
    }
    inline void EstimateBetweennessRanking(EstimateBetweenness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> EstimateBetweennessScores(EstimateBetweenness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<ForestCentrality> NewForestCentrality(
        const Graph &G, node root, double epsilon = 0.1, double kappa = 0.3)
    {
        return make_unique<ForestCentrality>(G, root, epsilon, kappa);
    }
    inline void ForestCentralityRanking(ForestCentrality &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ForestCentralityScores(ForestCentrality &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<double>> ForestCentralityGetDiagonal(const ForestCentrality &algo)
    {
        return make_unique<vector<double>>(algo.getDiagonal());
    }

    inline unique_ptr<GedWalk> NewGedWalk(
        const Graph &G, count k, double initEpsilon, double alpha,
        uint8_t bs, uint8_t gs,
        double spectralDelta)
    {
        GedWalk::BoundStrategy bs_;
        switch (bs)
        {
        case 0:
            bs_ = GedWalk::BoundStrategy::NO;
            break;
        case 1:
            bs_ = GedWalk::BoundStrategy::SPECTRAL;
            break;
        case 2:
            bs_ = GedWalk::BoundStrategy::GEOMETRIC;
            break;
        case 3:
            bs_ = GedWalk::BoundStrategy::ADAPTIVE_GEOMETRIC;
            break;
        default:
            break;
        }
        GedWalk::GreedyStrategy gs_;
        switch (gs)
        {
        case 0:
            gs_ = GedWalk::GreedyStrategy::LAZY;
            break;
        case 1:
            gs_ = GedWalk::GreedyStrategy::STOCHASTIC;
            break;
        default:
            break;
        }
        return make_unique<GedWalk>(G, k, initEpsilon, alpha, bs_, gs_, spectralDelta);
    }
    inline unique_ptr<vector<node>> GedWalkGroupMaxGedWalk(const GedWalk &algo)
    {
        return make_unique<vector<node>>(algo.groupMaxGedWalk());
    }
    inline double GedWalkScoreOfGroup(GedWalk &algo, rust::Slice<const node> group, double epsilon)
    {
        return algo.scoreOfGroup(group.begin(), group.end(), epsilon);
    }
}

#endif // NK_CENTRALITY_H
