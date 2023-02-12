
#ifndef NK_EMBEDDING_H
#define NK_EMBEDDING_H

#include "rust/cxx.h"
#include <networkit/embedding/Node2Vec.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<Node2Vec> NewNode2Vec(const Graph &G, double P = 1, double Q = 1, count L = 80, count N = 10, count D = 128)
    {
        return make_unique<Node2Vec>(G, P, Q, L, N, D);
    }

    inline void Node2VecGetFeatures(const Node2Vec &algo, rust::Vec<double> &ret)
    {
        for (auto &&row : algo.getFeatures())
        {
            for (auto &&el : row)
            {
                ret.push_back(el);
            }
        }
    }
}

#endif // NK_EMBEDDING_H