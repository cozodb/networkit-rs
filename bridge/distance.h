
#ifndef NK_DISTANCE_H
#define NK_DISTANCE_H

#include "rust/cxx.h"
#include "graph_event.h"
#include <networkit/distance/APSP.hpp>
#include <networkit/distance/AStar.hpp>
#include <networkit/distance/AdamicAdarDistance.hpp>
#include <networkit/distance/AlgebraicDistance.hpp>
#include <networkit/reachability/AllSimplePaths.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/BidirectionalBFS.hpp>
#include <networkit/distance/BidirectionalDijkstra.hpp>
#include <networkit/distance/CommuteTimeDistance.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/DynAPSP.hpp>
#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/Eccentricity.hpp>
#include <networkit/distance/EffectiveDiameter.hpp>
#include <networkit/distance/EffectiveDiameterApproximation.hpp>
#include <networkit/distance/HopPlotApproximation.hpp>
#include <networkit/distance/JaccardDistance.hpp>
#include <networkit/distance/MultiTargetBFS.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/distance/NeighborhoodFunction.hpp>
#include <networkit/distance/NeighborhoodFunctionApproximation.hpp>
#include <networkit/distance/NeighborhoodFunctionHeuristic.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/distance/ReverseBFS.hpp>
#include <networkit/distance/Volume.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<APSP> NewAPSP(const Graph &G)
    {
        return make_unique<APSP>(G);
    }

    inline node APSPGetDistances(const APSP &algo, rust::Vec<edgeweight> &ws)
    {
        auto distances = algo.getDistances();
        for (auto &&edges : distances)
        {
            for (auto &&wt : edges)
            {
                ws.push_back(wt);
            }
        }
        return distances.size();
    }

    inline unique_ptr<AStar> NewAStar(const Graph &G, const std::vector<double> &distanceHeu, node source, node target, bool storePred)
    {
        return make_unique<AStar>(G, distanceHeu, source, target, storePred);
    }

    inline unique_ptr<vector<node>> AStarGetPath(const AStar &algo)
    {
        return make_unique<vector<node>>(algo.getPath());
    }

    inline unique_ptr<vector<node>> AStarGetPredecessors(const AStar &algo)
    {
        return make_unique<vector<node>>(algo.getPredecessors());
    }

    inline unique_ptr<vector<edgeweight>> AStarGetDistances(const AStar &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }
    inline void AStarGetTargetIndexMap(const AStar &algo, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&pair : algo.getTargetIndexMap())
        {
            src.push_back(pair.first);
            dst.push_back(pair.second);
        }
    }
    inline void AStarSetTargets(AStar &algo, rust::Slice<const node> targets)
    {
        algo.setTargets(targets.begin(), targets.end());
    }

    inline unique_ptr<AdamicAdarDistance> NewAdamicAdarDistance(const Graph &G)
    {
        return make_unique<AdamicAdarDistance>(G);
    }

    inline unique_ptr<vector<double>> AdamicAdarDistanceGetEdgeScores(const AdamicAdarDistance &algo)
    {
        return make_unique<vector<double>>(algo.getEdgeScores());
    }

    inline unique_ptr<AlgebraicDistance> NewAlgebraicDistance(const Graph &G, count numberSystems = 10, count numberIterations = 30,
                                                              double omega = 0.5, index norm = 0, bool withEdgeScores = false)
    {
        return make_unique<AlgebraicDistance>(G, numberSystems, numberIterations, omega, norm, withEdgeScores);
    }

    inline unique_ptr<vector<double>> AlgebraicDistanceGetEdgeScores(const AlgebraicDistance &algo)
    {
        return make_unique<vector<double>>(algo.getEdgeScores());
    }

    inline unique_ptr<AllSimplePaths> NewAllSimplePaths(const Graph &G, node source, node target, count cutoff = none)
    {
        return make_unique<AllSimplePaths>(G, source, target, cutoff);
    }

    inline void AllSimplePathsGetAllSimplePaths(AllSimplePaths &algo, rust::Vec<node> &vs)
    {
        for (auto &&inner : algo.getAllSimplePaths())
        {
            for (auto &&j : inner)
            {
                vs.push_back(j);
            }
            vs.push_back(none);
        }
    }

    inline unique_ptr<BFS> NewBFS(const Graph &G, node source, bool storePaths = true, bool storeNodesSortedByDistance = false, node target = none)
    {
        return make_unique<BFS>(G, source, storePaths, storeNodesSortedByDistance, target);
    }

    inline unique_ptr<vector<edgeweight>> BFSGetDistances(BFS &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }

    inline unique_ptr<vector<node>> BFSGetPredecessors(const BFS &algo, node t)
    {
        return make_unique<vector<node>>(algo.getPredecessors(t));
    }

    inline unique_ptr<vector<node>> BFSGetPath(const BFS &algo, node t, bool forward)
    {
        return make_unique<vector<node>>(algo.getPath(t, forward));
    }

    inline void BFSGetPaths(const BFS &algo, node t, bool forward, rust::Vec<node> &vs)
    {
        for (auto &&p : algo.getPaths(t, forward))
        {
            vs.push_back(none);
            for (auto &&n : p)
            {
                vs.push_back(n);
            }
        }
    }

    inline unique_ptr<vector<node>> BFSGetNodeSortedByDistance(const BFS &algo)
    {
        return make_unique<vector<node>>(algo.getNodesSortedByDistance());
    }

    inline unique_ptr<BidirectionalBFS> NewBidirectionalBFS(const Graph &G, node source, node target, bool storePred)
    {
        return make_unique<BidirectionalBFS>(G, source, target, storePred);
    }

    inline unique_ptr<vector<node>> BidirectionalBFSGetPath(const BidirectionalBFS &algo)
    {
        return make_unique<vector<node>>(algo.getPath());
    }

    inline unique_ptr<vector<node>> BidirectionalBFSGetPredecessors(const BidirectionalBFS &algo)
    {
        return make_unique<vector<node>>(algo.getPredecessors());
    }

    inline unique_ptr<vector<edgeweight>> BidirectionalBFSGetDistances(const BidirectionalBFS &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }
    inline void BidirectionalBFSGetTargetIndexMap(const BidirectionalBFS &algo, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&pair : algo.getTargetIndexMap())
        {
            src.push_back(pair.first);
            dst.push_back(pair.second);
        }
    }
    inline void BidirectionalBFSSetTargets(BidirectionalBFS &algo, rust::Slice<const node> targets)
    {
        algo.setTargets(targets.begin(), targets.end());
    }

    inline unique_ptr<BidirectionalDijkstra> NewBidirectionalDijkstra(const Graph &G, node source, node target, bool storePred)
    {
        return make_unique<BidirectionalDijkstra>(G, source, target, storePred);
    }

    inline unique_ptr<vector<node>> BidirectionalDijkstraGetPath(const BidirectionalDijkstra &algo)
    {
        return make_unique<vector<node>>(algo.getPath());
    }

    inline unique_ptr<vector<node>> BidirectionalDijkstraGetPredecessors(const BidirectionalDijkstra &algo)
    {
        return make_unique<vector<node>>(algo.getPredecessors());
    }

    inline unique_ptr<vector<edgeweight>> BidirectionalDijkstraGetDistances(const BidirectionalDijkstra &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }
    inline void BidirectionalDijkstraGetTargetIndexMap(const BidirectionalDijkstra &algo, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&pair : algo.getTargetIndexMap())
        {
            src.push_back(pair.first);
            dst.push_back(pair.second);
        }
    }
    inline void BidirectionalDijkstraSetTargets(BidirectionalDijkstra &algo, rust::Slice<const node> targets)
    {
        algo.setTargets(targets.begin(), targets.end());
    }

    inline unique_ptr<CommuteTimeDistance> NewCommuteTimeDistance(const Graph &G, double tol = 0.1)
    {
        return make_unique<CommuteTimeDistance>(G, tol);
    }

    inline unique_ptr<Diameter> NewDiameter(const Graph &G, uint8_t algo_e, double error = -1.f, count nSamples = 0)
    {
        DiameterAlgo algo;
        switch (algo_e)
        {
        case 0:
            algo = DiameterAlgo::AUTOMATIC;
            break;
        case 1:
            algo = DiameterAlgo::EXACT;
            break;
        case 2:
            algo = DiameterAlgo::
                ESTIMATED_RANGE;
            break;
        case 3:
            algo = DiameterAlgo::
                ESTIMATED_SAMPLES;
            break;
        case 4:
            algo = DiameterAlgo::
                ESTIMATED_PEDANTIC;
            break;
        }
        return make_unique<Diameter>(G, algo, error, nSamples);
    }
    inline void DiameterGetDiameter(const Diameter &algo, count &lower, count &upper)
    {
        auto pair = algo.getDiameter();
        lower = pair.first;
        upper = pair.second;
    }

    inline unique_ptr<Dijkstra> NewDijkstra(const Graph &G, node source, bool storePaths = true, bool storeNodesSortedByDistance = false, node target = none)
    {
        return make_unique<Dijkstra>(G, source, storePaths, storeNodesSortedByDistance, target);
    }

    inline unique_ptr<vector<edgeweight>> DijkstraGetDistances(Dijkstra &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }

    inline unique_ptr<vector<node>> DijkstraGetPredecessors(const Dijkstra &algo, node t)
    {
        return make_unique<vector<node>>(algo.getPredecessors(t));
    }

    inline unique_ptr<vector<node>> DijkstraGetPath(const Dijkstra &algo, node t, bool forward)
    {
        return make_unique<vector<node>>(algo.getPath(t, forward));
    }

    inline void DijkstraGetPaths(const Dijkstra &algo, node t, bool forward, rust::Vec<node> &vs)
    {
        for (auto &&p : algo.getPaths(t, forward))
        {
            vs.push_back(none);
            for (auto &&n : p)
            {
                vs.push_back(n);
            }
        }
    }

    inline unique_ptr<vector<node>> DijkstraGetNodeSortedByDistance(const Dijkstra &algo)
    {
        return make_unique<vector<node>>(algo.getNodesSortedByDistance());
    }

    inline unique_ptr<DynAPSP> NewDynAPSP(Graph &G)
    {
        return make_unique<DynAPSP>(G);
    }

    inline node DynAPSPGetDistances(const DynAPSP &algo, rust::Vec<edgeweight> &ws)
    {
        auto distances = algo.getDistances();
        for (auto &&edges : distances)
        {
            for (auto &&wt : edges)
            {
                ws.push_back(wt);
            }
        }
        return distances.size();
    }

    inline void DynAPSPUpdate(DynAPSP &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynAPSPUpdateBatch(DynAPSP &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<DynBFS> NewDynBFS(const Graph &G, node s, bool storePredecessors = true)
    {
        return make_unique<DynBFS>(G, s, storePredecessors);
    }

    inline unique_ptr<vector<edgeweight>> DynBFSGetDistances(DynBFS &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }

    inline unique_ptr<vector<node>> DynBFSGetPredecessors(const DynBFS &algo, node t)
    {
        return make_unique<vector<node>>(algo.getPredecessors(t));
    }

    inline unique_ptr<vector<node>> DynBFSGetPath(const DynBFS &algo, node t, bool forward)
    {
        return make_unique<vector<node>>(algo.getPath(t, forward));
    }

    inline void DynBFSGetPaths(const DynBFS &algo, node t, bool forward, rust::Vec<node> &vs)
    {
        for (auto &&p : algo.getPaths(t, forward))
        {
            vs.push_back(none);
            for (auto &&n : p)
            {
                vs.push_back(n);
            }
        }
    }

    inline unique_ptr<vector<node>> DynBFSGetNodeSortedByDistance(const DynBFS &algo)
    {
        return make_unique<vector<node>>(algo.getNodesSortedByDistance());
    }

    inline void DynBFSUpdate(DynBFS &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynBFSUpdateBatch(DynBFS &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline void EccentricityGetValue(const Graph &G, node u, node &fartherest, count &dist)
    {
        auto res = Eccentricity::getValue(G, u);
        fartherest = res.first;
        dist = res.second;
    }

    inline unique_ptr<EffectiveDiameter> NewEffectiveDiameter(const Graph &G, double ratio)
    {
        return make_unique<EffectiveDiameter>(G, ratio);
    }

    inline unique_ptr<EffectiveDiameterApproximation> NewEffectiveDiameterApproximation(const Graph &G, double ratio, count k, count r)
    {
        return make_unique<EffectiveDiameterApproximation>(G, ratio, k, r);
    }

    inline unique_ptr<HopPlotApproximation> NewHopPlotApproximation(const Graph &G, count maxDistance, count k, count r)
    {
        return make_unique<HopPlotApproximation>(G, maxDistance, k, r);
    }

    inline void HopPlotApproximationGetHopPlot(const HopPlotApproximation &algo, rust::Vec<count> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.getHopPlot())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline unique_ptr<JaccardDistance> NewJaccardDistance(const Graph &G, const std::vector<count> &triangles)
    {
        return make_unique<JaccardDistance>(G, triangles);
    }

    inline unique_ptr<vector<double>> JaccardDistanceGetEdgeScores(const JaccardDistance &algo)
    {
        return make_unique<vector<double>>(algo.getEdgeScores());
    }

    inline unique_ptr<MultiTargetBFS> NewMultiTargetBFS(const Graph &G, node source, rust::Slice<const node> targets)
    {
        return make_unique<MultiTargetBFS>(G, source, targets.begin(), targets.end());
    }

    inline unique_ptr<vector<node>> MultiTargetBFSGetPath(const MultiTargetBFS &algo)
    {
        return make_unique<vector<node>>(algo.getPath());
    }

    inline unique_ptr<vector<node>> MultiTargetBFSGetPredecessors(const MultiTargetBFS &algo)
    {
        return make_unique<vector<node>>(algo.getPredecessors());
    }

    inline unique_ptr<vector<edgeweight>> MultiTargetBFSGetDistances(const MultiTargetBFS &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }
    inline void MultiTargetBFSGetTargetIndexMap(const MultiTargetBFS &algo, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&pair : algo.getTargetIndexMap())
        {
            src.push_back(pair.first);
            dst.push_back(pair.second);
        }
    }
    inline void MultiTargetBFSSetTargets(MultiTargetBFS &algo, rust::Slice<const node> targets)
    {
        algo.setTargets(targets.begin(), targets.end());
    }

    inline unique_ptr<MultiTargetDijkstra> NewMultiTargetDijkstra(const Graph &G, node source, rust::Slice<const node> targets)
    {
        return make_unique<MultiTargetDijkstra>(G, source, targets.begin(), targets.end());
    }

    inline unique_ptr<vector<node>> MultiTargetDijkstraGetPath(const MultiTargetDijkstra &algo)
    {
        return make_unique<vector<node>>(algo.getPath());
    }

    inline unique_ptr<vector<node>> MultiTargetDijkstraGetPredecessors(const MultiTargetDijkstra &algo)
    {
        return make_unique<vector<node>>(algo.getPredecessors());
    }

    inline unique_ptr<vector<edgeweight>> MultiTargetDijkstraGetDistances(const MultiTargetDijkstra &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }
    inline void MultiTargetDijkstraGetTargetIndexMap(const MultiTargetDijkstra &algo, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&pair : algo.getTargetIndexMap())
        {
            src.push_back(pair.first);
            dst.push_back(pair.second);
        }
    }
    inline void MultiTargetDijkstraSetTargets(MultiTargetDijkstra &algo, rust::Slice<const node> targets)
    {
        algo.setTargets(targets.begin(), targets.end());
    }

    inline unique_ptr<NeighborhoodFunction> NewNeighborhoodFunction(const Graph &G)
    {
        return make_unique<NeighborhoodFunction>(G);
    }

    inline unique_ptr<vector<count>> NeighborhoodFunctionGetNeighborhoodFunction(const NeighborhoodFunction &algo)
    {
        return make_unique<vector<count>>(algo.getNeighborhoodFunction());
    }

    inline unique_ptr<NeighborhoodFunctionApproximation> NewNeighborhoodFunctionApproximation(const Graph &G, count k = 64, count r = 7)
    {
        return make_unique<NeighborhoodFunctionApproximation>(G, k, r);
    }

    inline unique_ptr<vector<count>> NeighborhoodFunctionApproximationGetNeighborhoodFunction(const NeighborhoodFunctionApproximation &algo)
    {
        return make_unique<vector<count>>(algo.getNeighborhoodFunction());
    }

    inline unique_ptr<NeighborhoodFunctionHeuristic> NewNeighborhoodFunctionHeuristic(const Graph &G, count n_samples, uint8_t strategy)
    {
        NeighborhoodFunctionHeuristic::SelectionStrategy s;
        switch (strategy)
        {
        case 0:
            s = NeighborhoodFunctionHeuristic::SelectionStrategy::RANDOM;
            break;
        case 1:
            s = NeighborhoodFunctionHeuristic::SelectionStrategy::SPLIT;
            break;
        }
        return make_unique<NeighborhoodFunctionHeuristic>(G, n_samples, s);
    }

    inline unique_ptr<vector<count>> NeighborhoodFunctionHeuristicGetNeighborhoodFunction(const NeighborhoodFunctionHeuristic &algo)
    {
        return make_unique<vector<count>>(algo.getNeighborhoodFunction());
    }

    inline unique_ptr<PrunedLandmarkLabeling> NewPrunedLandmarkLabeling(const Graph &G)
    {
        return make_unique<PrunedLandmarkLabeling>(G);
    }

    inline unique_ptr<ReverseBFS> NewReverseBFS(const Graph &G, node source, bool storePaths = true, bool storeNodesSortedByDistance = false, node target = none)
    {
        return make_unique<ReverseBFS>(G, source, storePaths, storeNodesSortedByDistance, target);
    }

    inline unique_ptr<vector<edgeweight>> ReverseBFSGetDistances(ReverseBFS &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getDistances());
    }

    inline unique_ptr<vector<node>> ReverseBFSGetPredecessors(const ReverseBFS &algo, node t)
    {
        return make_unique<vector<node>>(algo.getPredecessors(t));
    }

    inline unique_ptr<vector<node>> ReverseBFSGetPath(const ReverseBFS &algo, node t, bool forward)
    {
        return make_unique<vector<node>>(algo.getPath(t, forward));
    }

    inline void ReverseBFSGetPaths(const ReverseBFS &algo, node t, bool forward, rust::Vec<node> &vs)
    {
        for (auto &&p : algo.getPaths(t, forward))
        {
            vs.push_back(none);
            for (auto &&n : p)
            {
                vs.push_back(n);
            }
        }
    }

    inline unique_ptr<vector<node>> ReverseBFSGetNodeSortedByDistance(const ReverseBFS &algo)
    {
        return make_unique<vector<node>>(algo.getNodesSortedByDistance());
    }

    inline unique_ptr<vector<double>> VolumeVolumes(const Graph &G, rust::Slice<const double> rs, count samples)
    {
        vector<double> radii{rs.begin(), rs.end()};
        return make_unique<vector<double>>(Volume::volume(G, radii, samples));
    }

    inline double VolumeVolume(const Graph &G, double r, count samples)
    {

        return Volume::volume(G, r, samples);
    }
}

#endif // NK_DISTANCE_H