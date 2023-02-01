//
// Created by Ziyang Hu on 2023/1/29.
//

#ifndef NETWORKIT_EXTRA_H
#define NETWORKIT_EXTRA_H

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/community/AdjustedRandMeasure.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include "rust/cxx.h"

namespace NetworKit
{

    using namespace std;

    // GRAPH

    inline unique_ptr<Graph> NewGraph(count n, bool weighted, bool directed, bool edgesIndexed)
    {
        return make_unique<Graph>(n, weighted, directed, edgesIndexed);
    }

    inline unique_ptr<Graph> CopyGraph(const Graph &g)
    {
        return make_unique<Graph>(g);
    }

    class GraphNodeIter
    {
        Graph::NodeIterator cur;
        Graph::NodeIterator end;

    public:
        GraphNodeIter(Graph::NodeIterator cur_, Graph::NodeIterator end_) : cur(cur_), end(end_) {}

        inline bool advance(node &u)
        {
            if (cur != end)
            {
                u = *cur;
                ++cur;
                return true;
            }
            else
            {
                return false;
            }
        }
    };

    inline unique_ptr<GraphNodeIter> NewGraphNodeIter(const Graph &g)
    {
        auto range = g.nodeRange();
        return make_unique<GraphNodeIter>(range.begin(), range.end());
    }

    class GraphEdgeIter
    {
        Graph::EdgeIterator cur;
        Graph::EdgeIterator end;

    public:
        GraphEdgeIter(Graph::EdgeIterator cur_, Graph::EdgeIterator end_) : cur(cur_), end(end_) {}

        inline bool advance(node &u, node &v)
        {
            if (cur != end)
            {
                u = (*cur).u;
                v = (*cur).v;
                ++cur;
                return true;
            }
            else
            {
                return false;
            }
        }
    };

    inline unique_ptr<GraphEdgeIter> NewGraphEdgeIter(const Graph &g)
    {
        auto range = g.edgeRange();
        return make_unique<GraphEdgeIter>(range.begin(), range.end());
    }

    class GraphEdgeWeightIter
    {
        Graph::EdgeWeightIterator cur;
        Graph::EdgeWeightIterator end;

    public:
        GraphEdgeWeightIter(Graph::EdgeWeightIterator cur_, Graph::EdgeWeightIterator end_) : cur(cur_), end(end_) {}

        inline bool advance(node &u, node &v, edgeweight &wt)
        {
            if (cur != end)
            {
                u = (*cur).u;
                v = (*cur).v;
                wt = (*cur).weight;
                ++cur;
                return true;
            }
            else
            {
                return false;
            }
        }
    };

    inline unique_ptr<GraphEdgeWeightIter> NewGraphEdgeWeightIter(const Graph &g)
    {
        auto range = g.edgeWeightRange();
        return make_unique<GraphEdgeWeightIter>(range.begin(), range.end());
    }

    class GraphNeighbourIter
    {
        Graph::NeighborIterator cur;
        Graph::NeighborIterator end;

    public:
        GraphNeighbourIter(Graph::NeighborIterator cur_, Graph::NeighborIterator end_) : cur(cur_), end(end_) {}

        inline bool advance(node &u)
        {
            if (cur != end)
            {
                u = *cur;
                ++cur;
                return true;
            }
            else
            {
                return false;
            }
        }
    };

    inline unique_ptr<GraphNeighbourIter> NewGraphNeighbourIter(const Graph &g, node u, bool in_neighbours)
    {
        if (in_neighbours)
        {
            auto range = g.inNeighborRange(u);
            return make_unique<GraphNeighbourIter>(range.begin(), range.end());
        }
        else
        {
            auto range = g.neighborRange(u);
            return make_unique<GraphNeighbourIter>(range.begin(), range.end());
        }
    }

    class GraphNeighbourWeightIter
    {
        Graph::NeighborWeightIterator cur;
        Graph::NeighborWeightIterator end;

    public:
        GraphNeighbourWeightIter(Graph::NeighborWeightIterator cur_, Graph::NeighborWeightIterator end_) : cur(cur_), end(end_) {}

        inline node current() const
        {
            return (*cur).first;
        }
        inline edgeweight current_weight() const
        {
            return (*cur).second;
        }
        inline bool advance(node &u, edgeweight &wt)
        {
            if (cur != end)
            {
                u = (*cur).first;
                wt = (*cur).second;
                ++cur;
                return true;
            }
            else
            {
                return false;
            }
        }
    };

    inline unique_ptr<GraphNeighbourWeightIter> NewGraphNeighbourWeightIter(const Graph &g, node u, bool in_neighbours)
    {
        if (in_neighbours)
        {
            auto range = g.weightInNeighborRange(u);
            return make_unique<GraphNeighbourWeightIter>(range.begin(), range.end());
        }
        else
        {
            auto range = g.weightNeighborRange(u);
            return make_unique<GraphNeighbourWeightIter>(range.begin(), range.end());
        }
    }

    // GRAPH BUILDER

    inline unique_ptr<GraphBuilder> NewGraphBuilder(count n, bool weighted, bool directed)
    {
        return make_unique<GraphBuilder>(n, weighted, directed);
    }

    inline unique_ptr<Graph> GraphBuilderCompleteGraph(GraphBuilder &builder, bool parallel)
    {
        return make_unique<Graph>(builder.completeGraph(parallel));
    }

    // GRAPH TOOLS

    namespace GraphTools
    {
        inline unique_ptr<Graph> GTCopyNodes(const Graph &G)
        {
            return make_unique<Graph>(copyNodes(G));
        }
        inline unique_ptr<Graph> GTCreateAugmentedGraph(const Graph &G, node &root)
        {
            auto ret = createAugmentedGraph(G);
            root = ret.second;
            return make_unique<Graph>(move(ret.first));
        }
        inline unique_ptr<Graph> GTGetCompactedGraph(const Graph &G, bool random)
        {
            unordered_map<node, node> map;
            if (random)
            {
                map = getRandomContinuousNodeIds(G);
            }
            else
            {
                map = getContinuousNodeIds(G);
            }
            return make_unique<Graph>(getCompactedGraph(G, map));
        }

        inline double GTInVolume(const Graph &G, rust::Slice<const node> nodes)
        {
            return inVolume(G, nodes.begin(), nodes.end());
        }

        inline void GTRandomEdge(const Graph &G, bool uniform, node &src, node &dst)
        {
            auto ret = randomEdge(G, uniform);
            src = ret.first;
            dst = ret.second;
        }
        inline void GTRandomEdges(const Graph &G, count n, rust::Vec<node> &src, rust::Vec<node> &dst)
        {
            auto ret = randomEdges(G, n);
            for (auto &&pair : ret)
            {
                src.push_back(pair.first);
                dst.push_back(pair.second);
            }
        }
        inline unique_ptr<vector<node>> GTRandomNodes(const Graph &G, count n)
        {
            return make_unique<vector<node>>(randomNodes(G, n));
        }
        inline void GTRemoveEdgesFromIsolatedSet(Graph &G, rust::Slice<const node> nodes)
        {
            removeEdgesFromIsolatedSet(G, nodes.begin(), nodes.end());
        }
        inline void GTSize(const Graph &G, count &n_nodes, count &n_edges)
        {
            auto ret = size(G);
            n_nodes = ret.first;
            n_edges = ret.second;
        }
        inline unique_ptr<Graph> GTSubgraphAndNeighborsFromNodes(const Graph &G, rust::Slice<const node> nodes, bool includeOutNeighbors = false, bool includeInNeighbors = false)
        {
            unordered_set<node> ns(nodes.begin(), nodes.end());
            return make_unique<Graph>(subgraphAndNeighborsFromNodes(G, ns, includeOutNeighbors, includeInNeighbors));
        }
        inline unique_ptr<Graph> GTSubgraphFromNodes(const Graph &G, rust::Slice<const node> nodes)
        {
            unordered_set<node> ns(nodes.begin(), nodes.end());
            return make_unique<Graph>(subgraphFromNodes(G, ns));
        }
        inline unique_ptr<Graph> GTToUndirected(const Graph &G)
        {
            return make_unique<Graph>(toUndirected(G));
        }
        inline unique_ptr<Graph> GTToUnweighted(const Graph &G)
        {
            return make_unique<Graph>(toUnweighted(G));
        }
        inline unique_ptr<Graph> GTToWeighted(const Graph &G)
        {
            return make_unique<Graph>(toWeighted(G));
        }
        inline unique_ptr<vector<node>> GTTopologicalSort(const Graph &G)
        {
            return make_unique<vector<node>>(topologicalSort(G));
        }
        inline unique_ptr<Graph> GTTranspose(const Graph &G)
        {
            return make_unique<Graph>(transpose(G));
        }
        inline double GTVolume(const Graph &G, rust::Slice<const node> nodes)
        {
            if (nodes.empty())
            {
                return volume(G);
            }
            else
            {
                return volume(G, nodes.begin(), nodes.end());
            }
        }
    }

    // PARTITION

    inline unique_ptr<Partition> NewPartition(index z)
    {
        return make_unique<Partition>(z);
    }

    inline unique_ptr<Partition> CopyPartition(const Partition &p)
    {
        return make_unique<Partition>(p);
    }

    inline unique_ptr<vector<count>> PTSubsetSizes(const Partition &p)
    {
        return make_unique<vector<count>>(p.subsetSizes());
    }

    inline void PTSubsetSizeMap(const Partition &p, rust::Vec<node> &ks, rust::Vec<node> &sz)
    {
        for (auto &&pair : p.subsetSizeMap())
        {
            ks.push_back(pair.first);
            sz.push_back(pair.second);
        }
    }
    inline void PTGetMembers(const Partition &p, index s, rust::Vec<node> &rs)
    {
        for (auto &&res : p.getMembers(s))
        {
            rs.push_back(res);
        }
    }
    inline void PTGetSubsetIds(const Partition &p, rust::Vec<node> &rs)
    {
        for (auto &&res : p.getSubsetIds())
        {
            rs.push_back(res);
        }
    }
    inline unique_ptr<string> PTGetName(const Partition &p)
    {
        return make_unique<string>(p.getName());
    }
    inline void PTSetName(Partition &p, rust::Str name)
    {
        string n(name);
        p.setName(n);
    }

    // COVER

    inline unique_ptr<Cover> NewCover()
    {
        return make_unique<Cover>();
    }

    inline unique_ptr<Cover> NewCoverWithSize(index z)
    {
        return make_unique<Cover>(z);
    }

    inline unique_ptr<Cover> NewCoverFromPartition(const Partition &p)
    {
        return make_unique<Cover>(p);
    }

    inline unique_ptr<Cover> CopyCover(const Cover &c)
    {
        return make_unique<Cover>(c);
    }

    inline void CVGetMembers(const Cover &p, index s, rust::Vec<node> &rs)
    {
        for (auto &&res : p.getMembers(s))
        {
            rs.push_back(res);
        }
    }
    inline void CVGetSubsetIds(const Cover &p, rust::Vec<node> &rs)
    {
        for (auto &&res : p.getSubsetIds())
        {
            rs.push_back(res);
        }
    }
    inline void CVSubsetSizeMap(const Cover &p, rust::Vec<node> &ks, rust::Vec<node> &sz)
    {
        for (auto &&pair : p.subsetSizeMap())
        {
            ks.push_back(pair.first);
            sz.push_back(pair.second);
        }
    }

    inline unique_ptr<vector<count>> CVSubsetSizes(const Cover &p)
    {
        return make_unique<vector<count>>(p.subsetSizes());
    }

    inline unique_ptr<vector<count>> CVSubsetsOf(const Cover &c, node e)
    {
        auto r = c.subsetsOf(e);
        return make_unique<vector<count>>(r.begin(), r.end());
    }

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
}

#endif // NETWORKIT_EXTRA_H
