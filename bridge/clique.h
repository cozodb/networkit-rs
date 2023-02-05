
#ifndef NK_CLIQUE_H
#define NK_CLIQUE_H

#include "rust/cxx.h"
#include <networkit/clique/MaximalCliques.hpp>

namespace NetworKit
{
    using namespace std;
    inline unique_ptr<MaximalCliques> NewMaximalCliques(const Graph &G, bool maximumOnly = false)
    {
        return make_unique<MaximalCliques>(G, maximumOnly);
    }

    inline void MaximalCliquesGetCliques(MaximalCliques &algo, rust::Vec<index> &cliques, rust::Vec<node> &ns)
    {
        index n = 0;
        for (auto &&clique : algo.getCliques())
        {
            for (auto &&node : clique)
            {
                cliques.push_back(n);
                ns.push_back(node);
            }
            n++;
        }
    }
}

#endif // NK_CLIQUE_H
