//
// Created by Ziyang Hu on 2023/1/29.
//

#ifndef NETWORKIT_EXTRA_H
#define NETWORKIT_EXTRA_H

#include <networkit/graph/Graph.hpp>
#include "rust/cxx.h"

namespace NetworKit
{

    using namespace std;

    inline unique_ptr<Graph> NewGraph(count n, bool weighted, bool directed, bool edgesIndexed)
    {
        return make_unique<Graph>(n, weighted, directed, edgesIndexed);
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

    unique_ptr<GraphNodeIter> NewGraphNodeIter(const Graph &g)
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

    unique_ptr<GraphEdgeIter> NewGraphEdgeIter(const Graph &g)
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

    unique_ptr<GraphEdgeWeightIter> NewGraphEdgeWeightIter(const Graph &g)
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

    unique_ptr<GraphNeighbourIter> NewGraphNeighbourIter(const Graph &g, node u, bool in_neighbours)
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

    unique_ptr<GraphNeighbourWeightIter> NewGraphNeighbourWeightIter(const Graph &g, node u, bool in_neighbours)
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

}

#endif // NETWORKIT_EXTRA_H
