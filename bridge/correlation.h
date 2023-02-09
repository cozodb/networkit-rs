
#ifndef NK_CORRELATION_H
#define NK_CORRELATION_H

#include "rust/cxx.h"
#include <networkit/correlation/Assortativity.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<Assortativity> NewAssortativity(const Graph &G, rust::Slice<const double> attrs)
    {
        vector<double> a{attrs.begin(), attrs.end()};
        return make_unique<Assortativity>(G, a);
    }

}

#endif // NK_CORRELATION_H