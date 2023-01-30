//
// Created by Ziyang Hu on 2023/1/29.
//

#ifndef NETWORKIT_EXTRA_H
#define NETWORKIT_EXTRA_H

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {
    PageRank NewPageRank(const Graph &G, double damp = 0.85, double tol = 1e-8, bool normalized = false,
                PageRank::SinkHandling distributeSinks = PageRank::SinkHandling::NO_SINK_HANDLING);
}

#endif //NETWORKIT_EXTRA_H
