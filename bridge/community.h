#ifndef NK_COMMUNITY_H
#define NK_COMMUNITY_H

#include <networkit/community/AdjustedRandMeasure.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/LocalCommunityEvaluation.hpp>
#include <networkit/community/CoverF1Similarity.hpp>
#include <networkit/community/CoverHubDominance.hpp>
#include <networkit/community/Coverage.hpp>
#include <networkit/community/CutClustering.hpp>
#include <networkit/community/EdgeCut.hpp>
#include <networkit/community/GraphClusteringTools.hpp>
#include <networkit/community/GraphStructuralRandMeasure.hpp>
#include <networkit/community/HubDominance.hpp>
#include <networkit/community/IntrapartitionDensity.hpp>
#include <networkit/community/IsolatedInterpartitionConductance.hpp>
#include <networkit/community/IsolatedInterpartitionExpansion.hpp>
#include <networkit/community/JaccardMeasure.hpp>
#include <networkit/community/LFM.hpp>
#include <networkit/community/LPDegreeOrdered.hpp>
#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/NMIDistance.hpp>
#include <networkit/community/NodeStructuralRandMeasure.hpp>
#include <networkit/community/OverlappingNMIDistance.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/community/ParallelLeiden.hpp>
#include <networkit/community/PartitionFragmentation.hpp>
#include <networkit/community/PartitionHubDominance.hpp>
#include <networkit/community/PartitionIntersection.hpp>
#include <networkit/community/StablePartitionNodes.hpp>
#include "rust/cxx.h"

namespace NetworKit
{

    using namespace std;

    // COMMUNITY

    inline unique_ptr<AdjustedRandMeasure> NewAdjustedRandMeasure()
    {
        return make_unique<AdjustedRandMeasure>();
    }

    inline unique_ptr<ClusteringGenerator> NewClusteringGenerator()
    {
        return make_unique<ClusteringGenerator>();
    }

    inline unique_ptr<Partition> CMMakeContinuousBalancedClustering(ClusteringGenerator &gen, const Graph &G, count k)
    {
        return make_unique<Partition>(gen.makeContinuousBalancedClustering(G, k));
    }

    inline unique_ptr<Partition> CMMakeNoncontinuousBalancedClustering(ClusteringGenerator &gen, const Graph &G, count k)
    {
        return make_unique<Partition>(gen.makeNoncontinuousBalancedClustering(G, k));
    }

    inline unique_ptr<Partition> CMMakeOneClustering(ClusteringGenerator &gen, const Graph &G)
    {
        return make_unique<Partition>(gen.makeOneClustering(G));
    }

    inline unique_ptr<Partition> CMMakeRandomClustering(ClusteringGenerator &gen, const Graph &G, count k)
    {
        return make_unique<Partition>(gen.makeRandomClustering(G, k));
    }

    inline unique_ptr<Partition> CMMakeSingletonClustering(ClusteringGenerator &gen, const Graph &G)
    {
        return make_unique<Partition>(gen.makeSingletonClustering(G));
    }

    inline unique_ptr<vector<double>> CoverF1SimilarityGetValues(const CoverF1Similarity &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<CoverF1Similarity> NewCoverF1Similarity(const Graph &G, const Cover &C, const Cover &reference)
    {
        return make_unique<CoverF1Similarity>(CoverF1Similarity(G, C, reference));
    }

    inline unique_ptr<vector<double>> CoverHubDominanceGetValues(const CoverHubDominance &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<CoverHubDominance> NewCoverHubDominance(const Graph &G, const Cover &C)
    {
        return make_unique<CoverHubDominance>(CoverHubDominance(G, C));
    }

    inline unique_ptr<Coverage> NewCoverage()
    {
        return make_unique<Coverage>();
    }

    inline unique_ptr<CutClustering> NewCutClustering(const Graph &G, edgeweight alpha)
    {
        return make_unique<CutClustering>(G, alpha);
    }

    inline unique_ptr<Partition> CutClusteringGetPartition(CutClustering &a)
    {
        return make_unique<Partition>(a.getPartition());
    }

    class HierarchyIter
    {
        map<edgeweight, Partition> data;
        map<edgeweight, Partition>::iterator cur;
        map<edgeweight, Partition>::iterator end;

    public:
        HierarchyIter(map<edgeweight, Partition> data_) : data(data_)
        {
            cur = data.begin();
            end = data.end();
        }
        inline bool isAtEnd() const
        {
            return cur == end;
        }
        inline void advance()
        {
            if (cur != end)
            {
                cur++;
            }
        }
        inline edgeweight curKey() const
        {
            return cur->first;
        }
        inline unique_ptr<Partition> curVal() const
        {
            return make_unique<Partition>(cur->second);
        }
    };

    // CutClustering::getClusterHierarchy
    inline unique_ptr<HierarchyIter> CutClusteringGetClusterHierarchy(const Graph &g)
    {
        return make_unique<HierarchyIter>(CutClustering::getClusterHierarchy(g));
    }

    inline unique_ptr<EdgeCut> NewEdgeCut()
    {
        return make_unique<EdgeCut>();
    }

    namespace GraphClusteringTools
    {
        inline unique_ptr<Graph> MakeCommunicationGraph(const Graph &graph, Partition &zeta)
        {
            return make_unique<Graph>(communicationGraph(graph, zeta));
        }
    }

    inline unique_ptr<GraphStructuralRandMeasure> NewGraphStructuralRandMeasure()
    {
        return make_unique<GraphStructuralRandMeasure>();
    }

    inline unique_ptr<HubDominance> NewHubDominance()
    {
        return make_unique<HubDominance>();
    }

    inline unique_ptr<IntrapartitionDensity> NewIntrapartitionDensity(const Graph &G, const Partition &P)
    {
        return make_unique<IntrapartitionDensity>(G, P);
    }

    inline unique_ptr<vector<double>> IntrapartitionDensityGetValues(const IntrapartitionDensity &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<IsolatedInterpartitionConductance> NewIsolatedInterpartitionConductance(const Graph &G, const Partition &P)
    {
        return make_unique<IsolatedInterpartitionConductance>(G, P);
    }

    inline unique_ptr<vector<double>> IsolatedInterpartitionConductanceGetValues(const IsolatedInterpartitionConductance &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<IsolatedInterpartitionExpansion> NewIsolatedInterpartitionExpansion(const Graph &G, const Partition &P)
    {
        return make_unique<IsolatedInterpartitionExpansion>(G, P);
    }

    inline unique_ptr<vector<double>> IsolatedInterpartitionExpansionGetValues(const IsolatedInterpartitionExpansion &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<JaccardMeasure> NewJaccardMeasure()
    {
        return make_unique<JaccardMeasure>();
    }

    inline unique_ptr<LFM> NewLFM(const Graph &G, SelectiveCommunityDetector &scd)
    {
        return make_unique<LFM>(G, scd);
    }

    inline unique_ptr<Cover> LFMGetCover(const LFM &algo)
    {
        return make_unique<Cover>(algo.getCover());
    }

    inline unique_ptr<LPDegreeOrdered> NewLPDegreeOrdered(const Graph &g)
    {
        return make_unique<LPDegreeOrdered>(g);
    }

    inline unique_ptr<Partition> LPDegreeOrderedGetPartition(LPDegreeOrdered &a)
    {
        return make_unique<Partition>(a.getPartition());
    }

    inline unique_ptr<LouvainMapEquation> NewLouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations, const rust::Str parallelizationStrategy)
    {
        string ps(parallelizationStrategy);
        return make_unique<LouvainMapEquation>(graph, hierarchical, maxIterations, ps);
    }

    inline unique_ptr<Partition> LouvainMapEquationGetPartition(LouvainMapEquation &a)
    {
        return make_unique<Partition>(a.getPartition());
    }

    inline unique_ptr<Modularity> NewModularity()
    {
        return make_unique<Modularity>();
    }

    inline unique_ptr<NMIDistance> NewNMIDistance()
    {
        return make_unique<NMIDistance>();
    }

    inline unique_ptr<NodeStructuralRandMeasure> NewNodeStructuralRandMeasure()
    {
        return make_unique<NodeStructuralRandMeasure>();
    }

    inline unique_ptr<OverlappingNMIDistance> NewOverlappingNMIDistance(
        uint8_t normalization)
    {
        OverlappingNMIDistance::Normalization n;
        switch (normalization)
        {
        case 2:
            n = OverlappingNMIDistance::Normalization::ARITHMETIC_MEAN;
            break;
        case 4:
            n = OverlappingNMIDistance::Normalization::JOINT_ENTROPY;
            break;
        case 1:
            n = OverlappingNMIDistance::Normalization::GEOMETRIC_MEAN;
            break;
        case 3:
            n = OverlappingNMIDistance::Normalization::MAX;
            break;
        case 0:
            n = OverlappingNMIDistance::Normalization::MIN;
            break;
        default:
            n = OverlappingNMIDistance::Normalization::MIN;
            break;
        }
        return make_unique<OverlappingNMIDistance>(n);
    }

    inline unique_ptr<PLM> NewPLM(const Graph &G, bool refine, double gamma, rust::Str par,
                                  count maxIter, bool turbo, bool recurse)
    {
        return make_unique<PLM>(G, refine, gamma, string(par), maxIter, turbo, recurse);
    }

    inline unique_ptr<Graph> PLMCoarsen(const Graph &G, const Partition &zeta, rust::Vec<node> &mapping)
    {
        auto res = PLM::coarsen(G, zeta);
        for (auto &&n : res.second)
        {
            mapping.push_back(n);
        }
        return make_unique<Graph>(res.first);
    }
    inline unique_ptr<Partition> PLMProlong(
        const Graph &g, const Partition &zetaCoarse, const Graph &Gfine,
        rust::Slice<const node> nodeToMetaNode)
    {
        return make_unique<Partition>(PLM::prolong(g, zetaCoarse, Gfine, vector<node>(nodeToMetaNode.begin(), nodeToMetaNode.end())));
    }

    inline unique_ptr<Partition> PLMGetPartition(PLM &a)
    {
        return make_unique<Partition>(a.getPartition());
    }

    inline unique_ptr<PLP> NewPLP(const Graph &G, count theta = none, count maxIterations = none)
    {
        return make_unique<PLP>(G, theta, maxIterations);
    }

    inline unique_ptr<Partition> PLPGetPartition(PLP &a)
    {
        return make_unique<Partition>(a.getPartition());
    }

    inline unique_ptr<ParallelLeiden> NewParallelLeiden(const Graph &graph, count iterations = 3, bool randomize = true,
                                                        double gamma = 1)
    {
        return make_unique<ParallelLeiden>(graph, iterations, randomize, gamma);
    }

    inline unique_ptr<Partition> ParallelLeidenGetPartition(ParallelLeiden &a)
    {
        return make_unique<Partition>(a.getPartition());
    }

    inline unique_ptr<PartitionFragmentation> NewPartitionFragmentation(const Graph &G, const Partition &P)
    {
        return make_unique<PartitionFragmentation>(G, P);
    }

    inline unique_ptr<vector<double>> PartitionFragmentationGetValues(const PartitionFragmentation &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<PartitionHubDominance> NewPartitionHubDominance(const Graph &G, const Partition &P)
    {
        return make_unique<PartitionHubDominance>(G, P);
    }

    inline unique_ptr<vector<double>> PartitionHubDominanceGetValues(const PartitionHubDominance &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

    inline unique_ptr<PartitionIntersection> NewPartitionIntersection()
    {
        return make_unique<PartitionIntersection>();
    }

    inline unique_ptr<Partition> PartitionIntersectionCalculate(PartitionIntersection &algo, const Partition &zeta, const Partition &eta)
    {
        return make_unique<Partition>(algo.calculate(zeta, eta));
    }

    inline unique_ptr<StablePartitionNodes> NewStablePartitionNodes(const Graph &G, const Partition &P)
    {
        return make_unique<StablePartitionNodes>(G, P);
    }

    inline unique_ptr<vector<double>> StablePartitionNodesGetValues(const StablePartitionNodes &e)
    {
        return make_unique<vector<double>>(e.getValues());
    }

}

#endif // NK_COMMUNITY_H