
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

}

#endif // NK_CENTRALITY_H
