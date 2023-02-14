
#ifndef NK_MISC_H
#define NK_MISC_H

#include "rust/cxx.h"
#include <networkit/embedding/Node2Vec.hpp>
#include <networkit/flow/EdmondsKarp.hpp>
#include <networkit/independentset/Luby.hpp>

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

    inline unique_ptr<Luby> NewLuby()
    {
        return make_unique<Luby>();
    }

    inline void LubyRun(Luby &algo, const Graph &g, rust::Vec<bool> &ret)
    {
        for (auto &&v : algo.run(g))
        {
            ret.push_back(v);
        }
    }
    inline bool LubyIsIndependentSet(const Luby &algo, rust::Slice<const bool> set, const Graph &g)
    {
        vector<bool> v{set.begin(), set.end()};
        return algo.isIndependentSet(v, g);
    }
}

#endif // NK_MISC_H