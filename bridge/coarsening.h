#ifndef NK_COARSENING_H
#define NK_COARSENING_H

#include "rust/cxx.h"
#include <networkit/coarsening/GraphCoarsening.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<ParallelPartitionCoarsening> NewParallelPartitionCoarsening(const Graph &G, const Partition &zeta, bool parallel = true)
    {
        return make_unique<ParallelPartitionCoarsening>(G, zeta, parallel);
    }

    inline unique_ptr<Graph> ParallelPartitionCoarseningGetCoarseGraph(const ParallelPartitionCoarsening &algo)
    {
        return make_unique<Graph>(algo.getCoarseGraph());
    }

    inline unique_ptr<vector<node>> ParallelPartitionCoarseningGetFineToCoarseNodeMapping(const ParallelPartitionCoarsening &algo)
    {
        return make_unique<vector<node>>(algo.getFineToCoarseNodeMapping());
    }
}

#endif // NK_COARSENING_H