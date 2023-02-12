
#ifndef NK_GENERATORS_H
#define NK_GENERATORS_H

#include "rust/cxx.h"
#include <networkit/generators/BarabasiAlbertGenerator.hpp>

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
}

#endif // NK_GENERATORS_H