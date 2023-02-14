
#ifndef NK_GENERATORS_H
#define NK_GENERATORS_H

#include "rust/cxx.h"
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/generators/ChungLuGenerator.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicForestFireGenerator.hpp>
#include <networkit/generators/DynamicHyperbolicGenerator.hpp>
#include <networkit/generators/DynamicPathGenerator.hpp>
#include <networkit/generators/DynamicPubWebGenerator.hpp>
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/LFRGenerator.hpp>
#include <networkit/generators/MocnikGenerator.hpp>
#include <networkit/generators/MocnikGeneratorBasic.hpp>
#include <networkit/generators/PowerlawDegreeSequence.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/RmatGenerator.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<BarabasiAlbertGenerator> NewBarabasiAlbertGenerator(count k, count nMax, count n0 = 0, bool batagelj = true)
    {
        return make_unique<BarabasiAlbertGenerator>(k, nMax, n0, batagelj);
    }

    inline unique_ptr<Graph> BarabasiAlbertGeneratorGenerate(BarabasiAlbertGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<ChungLuGenerator> NewChungLuGenerator(const std::vector<count> &degreeSequence)
    {
        return make_unique<ChungLuGenerator>(degreeSequence);
    }

    inline unique_ptr<Graph> ChungLuGeneratorGenerate(ChungLuGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<ClusteredRandomGraphGenerator> NewClusteredRandomGraphGenerator(count n, count k, double pIntra, double pInter)
    {
        return make_unique<ClusteredRandomGraphGenerator>(n, k, pIntra, pInter);
    }

    inline unique_ptr<Graph> ClusteredRandomGraphGeneratorGenerate(ClusteredRandomGraphGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<Partition> ClusteredRandomGraphGeneratorGetCommunities(ClusteredRandomGraphGenerator &algo)
    {
        return make_unique<Partition>(algo.getCommunities());
    }

    inline unique_ptr<DorogovtsevMendesGenerator> NewDorogovtsevMendesGenerator(count nNodes)
    {
        return make_unique<DorogovtsevMendesGenerator>(nNodes);
    }

    inline unique_ptr<Graph> DorogovtsevMendesGeneratorGenerate(DorogovtsevMendesGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline uint8_t convertEventType(GraphEvent::Type t)
    {
        switch (t)
        {
        case GraphEvent::Type::NODE_ADDITION:
            return 0;
            break;
        case GraphEvent::Type::NODE_REMOVAL:
            return 1;
            break;
        case GraphEvent::Type::NODE_RESTORATION:
            return 2;
            break;
        case GraphEvent::Type::EDGE_ADDITION:
            return 3;
            break;
        case GraphEvent::Type::EDGE_REMOVAL:
            return 4;
            break;
        case GraphEvent::Type::EDGE_WEIGHT_UPDATE:
            return 5;
            break;
        case GraphEvent::Type::EDGE_WEIGHT_INCREMENT:
            return 6;
            break;
        case GraphEvent::Type::TIME_STEP:
            return 7;
            break;

        default:
            return 255;
            break;
        }
    }

    inline unique_ptr<DynamicDorogovtsevMendesGenerator> NewDynamicDorogovtsevMendesGenerator()
    {
        return make_unique<DynamicDorogovtsevMendesGenerator>();
    }

    inline void DynamicDorogovtsevMendesGeneratorGenerate(DynamicDorogovtsevMendesGenerator &algo, count nSteps, rust::Vec<uint8_t> &tps, rust::Vec<node> &us, rust::Vec<node> &vs, rust::Vec<edgeweight> &ws)
    {
        for (auto &&evt : algo.generate(nSteps))
        {
            tps.push_back(convertEventType(evt.type));
            us.push_back(evt.u);
            vs.push_back(evt.v);
            ws.push_back(evt.w);
        }
    }

    inline unique_ptr<DynamicForestFireGenerator> NewDynamicForestFireGenerator(double p, bool directed, double r = 1.0)
    {
        return make_unique<DynamicForestFireGenerator>(p, directed, r);
    }

    inline void DynamicForestFireGeneratorGenerate(DynamicForestFireGenerator &algo, count nSteps, rust::Vec<uint8_t> &tps, rust::Vec<node> &us, rust::Vec<node> &vs, rust::Vec<edgeweight> &ws)
    {
        for (auto &&evt : algo.generate(nSteps))
        {
            tps.push_back(convertEventType(evt.type));
            us.push_back(evt.u);
            vs.push_back(evt.v);
            ws.push_back(evt.w);
        }
    }

    inline unique_ptr<DynamicHyperbolicGenerator> NewDynamicHyperbolicGenerator(count n = 1000, double avgDegree = 6, double exp = 3, double T = 0,
                                                                                double moveEachStep = 0, double moveDistance = 0)
    {
        return make_unique<DynamicHyperbolicGenerator>(n, avgDegree, exp, T, moveEachStep, moveDistance);
    }

    inline void DynamicHyperbolicGeneratorGenerate(DynamicHyperbolicGenerator &algo, count nSteps, rust::Vec<uint8_t> &tps, rust::Vec<node> &us, rust::Vec<node> &vs, rust::Vec<edgeweight> &ws)
    {
        for (auto &&evt : algo.generate(nSteps))
        {
            tps.push_back(convertEventType(evt.type));
            us.push_back(evt.u);
            vs.push_back(evt.v);
            ws.push_back(evt.w);
        }
    }

    inline unique_ptr<Graph> DynamicHyperbolicGeneratorGetGraph(const DynamicHyperbolicGenerator &algo)
    {
        return make_unique<Graph>(algo.getGraph());
    }

    inline void DynamicHyperbolicGeneratorGetCoordinates(const DynamicHyperbolicGenerator &algo, rust::Vec<double> &xs, rust::Vec<double> &ys)
    {
        for (auto &&pair : algo.getCoordinates())
        {
            xs.push_back(pair[0]);
            ys.push_back(pair[1]);
        }
    }

    inline unique_ptr<DynamicPathGenerator> NewDynamicPathGenerator()
    {
        return make_unique<DynamicPathGenerator>();
    }

    inline void DynamicPathGeneratorGenerate(DynamicPathGenerator &algo, count nSteps, rust::Vec<uint8_t> &tps, rust::Vec<node> &us, rust::Vec<node> &vs, rust::Vec<edgeweight> &ws)
    {
        for (auto &&evt : algo.generate(nSteps))
        {
            tps.push_back(convertEventType(evt.type));
            us.push_back(evt.u);
            vs.push_back(evt.v);
            ws.push_back(evt.w);
        }
    }

    inline unique_ptr<DynamicPubWebGenerator> NewDynamicPubWebGenerator(count numNodes, count numberOfDenseAreas, coordinate neighborhoodRadius,
                                                                        count maxNumberOfNeighbors, bool writeInitialGraphToStream = true)
    {
        return make_unique<DynamicPubWebGenerator>(numNodes, numberOfDenseAreas, neighborhoodRadius,
                                                   maxNumberOfNeighbors, writeInitialGraphToStream);
    }

    inline void DynamicPubWebGeneratorGenerate(DynamicPubWebGenerator &algo, count nSteps, rust::Vec<uint8_t> &tps, rust::Vec<node> &us, rust::Vec<node> &vs, rust::Vec<edgeweight> &ws)
    {
        for (auto &&evt : algo.generate(nSteps))
        {
            tps.push_back(convertEventType(evt.type));
            us.push_back(evt.u);
            vs.push_back(evt.v);
            ws.push_back(evt.w);
        }
    }

    inline unique_ptr<Graph> DynamicPubWebGeneratorGetGraph(const DynamicPubWebGenerator &algo)
    {
        return make_unique<Graph>(algo.getGraph());
    }

    inline void DynamicPubWebGeneratorGetCoordinates(const DynamicPubWebGenerator &algo, rust::Vec<double> &xs, rust::Vec<double> &ys)
    {
        for (auto &&pair : algo.getCoordinates())
        {
            xs.push_back(pair[0]);
            ys.push_back(pair[1]);
        }
    }

    inline void DynamicPubWebGeneratorGetNewCoordinates(const DynamicPubWebGenerator &algo, rust::Vec<node> &ns, rust::Vec<double> &xs, rust::Vec<double> &ys)
    {
        for (auto &&pair : algo.getNewCoordinates())
        {
            ns.push_back(pair.first);
            xs.push_back(pair.second[0]);
            ys.push_back(pair.second[1]);
        }
    }

    inline unique_ptr<EdgeSwitchingMarkovChainGenerator> NewEdgeSwitchingMarkovChainGenerator(const std::vector<count> &sequence, bool ignoreIfNotRealizable = false, count numSwitchesPerEdge = 10)
    {
        return make_unique<EdgeSwitchingMarkovChainGenerator>(sequence, ignoreIfNotRealizable, numSwitchesPerEdge);
    }

    inline unique_ptr<Graph> EdgeSwitchingMarkovChainGeneratorGenerate(EdgeSwitchingMarkovChainGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<ErdosRenyiGenerator> NewErdosRenyiGenerator(count nNodes, double prob, bool directed = false, bool self_loops = false)
    {
        return make_unique<ErdosRenyiGenerator>(nNodes, prob, directed, self_loops);
    }

    inline unique_ptr<Graph> ErdosRenyiGeneratorGenerate(ErdosRenyiGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<HavelHakimiGenerator> NewHavelHakimiGenerator(const std::vector<count> &sequence, bool ignoreIfRealizable = false)
    {
        return make_unique<HavelHakimiGenerator>(sequence, ignoreIfRealizable);
    }

    inline unique_ptr<Graph> HavelHakimiGeneratorGenerate(HavelHakimiGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<HyperbolicGenerator> NewHyperbolicGenerator(count n = 10000, double avgDegree = 6, double exp = 3, double T = 0)
    {
        return make_unique<HyperbolicGenerator>(n, avgDegree, exp, T);
    }

    inline unique_ptr<Graph> HyperbolicGeneratorGenerate(HyperbolicGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<Graph> HyperbolicGeneratorGenerate(HyperbolicGenerator &algo, rust::Slice<const double> angles, rust::Slice<const double> radii, double R,
                                                         double T = 0)
    {
        vector<double> as{angles.begin(), angles.end()};
        vector<double> rs{radii.begin(), radii.end()};
        return make_unique<Graph>(algo.generate(as, rs, R, T));
    }

    inline unique_ptr<LFRGenerator> NewLFRGenerator(count n)
    {
        return make_unique<LFRGenerator>(n);
    }

    inline unique_ptr<Graph> LFRGeneratorGenerate(LFRGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline void LFRGeneratorSetDegreeSequence(LFRGenerator &algo, rust::Slice<const count> seq)
    {
        vector<count> s{seq.begin(), seq.end()};
        algo.setDegreeSequence(s);
    }
    inline void LFRGeneratorSetCommunitySizeSequence(LFRGenerator &algo, rust::Slice<const count> seq)
    {
        vector<count> s{seq.begin(), seq.end()};
        algo.setCommunitySizeSequence(s);
    }
    inline void LFRGeneratorSetPartition(LFRGenerator &algo, unique_ptr<Partition> p)
    {
        algo.setPartition(*p);
    }
    inline unique_ptr<Graph> LFRGeneratorGetGraph(const LFRGenerator &algo)
    {
        return make_unique<Graph>(algo.getGraph());
    }
    inline unique_ptr<Partition> LFRGeneratorGetPartition(const LFRGenerator &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline unique_ptr<MocnikGenerator> NewMocnikGenerator(count dim, count n, double k, bool weighted = false)
    {
        return make_unique<MocnikGenerator>(dim, n, k, weighted);
    }

    inline unique_ptr<Graph> MocnikGeneratorGenerate(MocnikGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<MocnikGeneratorBasic> NewMocnikGeneratorBasic(count dim, count n, double k)
    {
        return make_unique<MocnikGeneratorBasic>(dim, n, k);
    }

    inline unique_ptr<Graph> MocnikGeneratorBasicGenerate(MocnikGeneratorBasic &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<PowerlawDegreeSequence> NewPowerlawDegreeSequence(count minDeg, count maxDeg, double gamma)
    {
        return make_unique<PowerlawDegreeSequence>(minDeg, maxDeg, gamma);
    }
    inline unique_ptr<vector<count>> PowerlawDegreeSequenceGetDegreeSequence(const PowerlawDegreeSequence &algo, count numNodes)
    {
        return make_unique<vector<count>>(algo.getDegreeSequence(numNodes));
    }

    inline unique_ptr<PubWebGenerator> NewPubWebGenerator(count numNodes, count numberOfDenseAreas, coordinate neighborhoodRadius,
                                                          count maxNumberOfNeighbors)
    {
        return make_unique<PubWebGenerator>(numNodes, numberOfDenseAreas, neighborhoodRadius,
                                            maxNumberOfNeighbors);
    }

    inline unique_ptr<Graph> PubWebGeneratorGenerate(PubWebGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline void PubWebGeneratorGetCoordinates(const PubWebGenerator &algo, rust::Vec<double> &xs, rust::Vec<double> &ys)
    {
        for (auto &&pair : algo.getCoordinates())
        {
            xs.push_back(pair[0]);
            ys.push_back(pair[1]);
        }
    }

    inline unique_ptr<RegularRingLatticeGenerator> NewRegularRingLatticeGenerator(count nNodes, count nNeighbors)
    {
        return make_unique<RegularRingLatticeGenerator>(nNodes, nNeighbors);
    }

    inline unique_ptr<Graph> RegularRingLatticeGeneratorGenerate(RegularRingLatticeGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<RmatGenerator> NewRmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d,
                                                      bool weighted = false, count reduceNodes = 0)
    {
        return make_unique<RmatGenerator>(scale, edgeFactor, a, b, c, d, weighted, reduceNodes);
    }

    inline unique_ptr<Graph> RmatGeneratorGenerate(RmatGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<WattsStrogatzGenerator> NewWattsStrogatzGenerator(count nNodes, count nNeighbors, double p)
    {
        return make_unique<WattsStrogatzGenerator>(nNodes, nNeighbors, p);
    }

    inline unique_ptr<Graph> WattsStrogatzGeneratorGenerate(WattsStrogatzGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }
}

#endif // NK_GENERATORS_H