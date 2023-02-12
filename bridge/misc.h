
#ifndef NK_MISC_H
#define NK_MISC_H

#include "rust/cxx.h"
#include <networkit/embedding/Node2Vec.hpp>
#include <networkit/flow/EdmondsKarp.hpp>

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

    inline unique_ptr<EdmondsKarp> NewEdmondsKarp(const Graph &G, node source, node sink)
    {
        return make_unique<EdmondsKarp>(G, source, sink);
    }

    inline unique_ptr<vector<node>> EdmondsKarpGetSourceSet(const EdmondsKarp &algo)
    {
        return make_unique<vector<node>>(algo.getSourceSet());
    }

    inline unique_ptr<vector<edgeweight>> EdmondsKarpGetFlowVector(const EdmondsKarp &algo)
    {
        return make_unique<vector<edgeweight>>(algo.getFlowVector());
    }
}

#endif // NK_MISC_H