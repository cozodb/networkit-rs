//
// Created by Ziyang Hu on 2023/1/29.
//

#ifndef NETWORKIT_EXTRA_H
#define NETWORKIT_EXTRA_H

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

    using namespace std;

    inline unique_ptr <Graph> NewGraph(count n, bool weighted, bool directed, bool edgesIndexed) {
        return make_unique<Graph>(n, weighted, directed, edgesIndexed);
    }
//    using Graph_EdgeRange = Graph::EdgeRange;
//
//    inline unique_ptr<Graph_EdgeRange> NewEdgeRange(const Graph& g) {
//        return make_unique<Graph_EdgeRange>(g.edgeRange());
//    }

//    using Graph_NodeRange = Graph::NodeRange;
//    using Graph_NodeIterator = Graph::NodeIterator;

//    inline unique_ptr<Graph_NodeRange> NewNodeRange(const Graph& g) {
//        return make_unique<Graph_NodeRange>(g.nodeRange());
//    }
//    inline unique_ptr<
}

#endif //NETWORKIT_EXTRA_H
