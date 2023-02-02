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
        unique_ptr<Graph> MakeCommunicationGraph(const Graph &graph, Partition &zeta)
        {
            return make_unique<Graph>(communicationGraph(graph, zeta));
        }
    }
}

#endif // NK_COMMUNITY_H