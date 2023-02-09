
#ifndef NK_COMPONENT_H
#define NK_COMPONENT_H

#include "rust/cxx.h"
#include "graph_event.h"
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/DynConnectedComponents.hpp>
#include <networkit/components/DynWeaklyConnectedComponents.hpp>
#include <networkit/components/ParallelConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/components/WeaklyConnectedComponents.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<BiconnectedComponents> NewBiconnectedComponents(const Graph &G)
    {
        return make_unique<BiconnectedComponents>(G);
    }

    inline void BiconnectedComponentsGetComponentSizes(const BiconnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void BiconnectedComponentsGetComponents(const BiconnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }

    inline void BiconnectedComponentsGetComponentOfNode(const BiconnectedComponents &algo, node u, rust::Vec<node> &vs)
    {
        for (auto &&n : algo.getComponentsOfNode(u))
        {
            vs.push_back(n);
        }
    }

    inline unique_ptr<ConnectedComponents> NewConnectedComponents(const Graph &G)
    {
        return make_unique<ConnectedComponents>(G);
    }

    inline unique_ptr<Graph> ConnectedComponentsExtractLargestConnectedComponent(const Graph &G, bool compactGraph = false)
    {
        return make_unique<Graph>(ConnectedComponents::extractLargestConnectedComponent(G, compactGraph));
    }

    inline unique_ptr<Partition> ConnectedComponentsGetPartition(const ConnectedComponents &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline void ConnectedComponentsGetComponentSizes(const ConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void ConnectedComponentsGetComponents(const ConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }

    inline unique_ptr<DynConnectedComponents> NewDynConnectedComponents(const Graph &G)
    {
        return make_unique<DynConnectedComponents>(G);
    }

    inline unique_ptr<Partition> DynConnectedComponentsGetPartition(const DynConnectedComponents &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline void DynConnectedComponentsGetComponentSizes(const DynConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void DynConnectedComponentsGetComponents(const DynConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }

    inline void DynConnectedComponentsUpdate(DynConnectedComponents &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynConnectedComponentsUpdateBatch(DynConnectedComponents &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<DynWeaklyConnectedComponents> NewDynWeaklyConnectedComponents(const Graph &G)
    {
        return make_unique<DynWeaklyConnectedComponents>(G);
    }

    inline unique_ptr<Partition> DynWeaklyConnectedComponentsGetPartition(const DynWeaklyConnectedComponents &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline void DynWeaklyConnectedComponentsGetComponentSizes(const DynWeaklyConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void DynWeaklyConnectedComponentsGetComponents(const DynWeaklyConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }

    inline void DynWeaklyConnectedComponentsUpdate(DynWeaklyConnectedComponents &algo, uint8_t kind, node u, node v, edgeweight ew)
    {
        algo.update(toGraphEvent(kind, u, v, ew));
    }

    inline void DynWeaklyConnectedComponentsUpdateBatch(DynWeaklyConnectedComponents &algo, rust::Slice<const uint8_t> kinds, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const edgeweight> ews)
    {
        vector<GraphEvent> evs;
        evs.reserve(kinds.length());
        for (size_t i = 0; i < kinds.length(); ++i)
        {
            evs.emplace_back(toGraphEvent(kinds[i], us[i], vs[i], ews[i]));
        }
        algo.updateBatch(evs);
    }

    inline unique_ptr<ParallelConnectedComponents> NewParallelConnectedComponents(const Graph &G, bool coarsening = true)
    {
        return make_unique<ParallelConnectedComponents>(G, coarsening);
    }

    inline unique_ptr<Partition> ParallelConnectedComponentsGetPartition(const ParallelConnectedComponents &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline void ParallelConnectedComponentsGetComponentSizes(const ParallelConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void ParallelConnectedComponentsGetComponents(const ParallelConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }

    inline unique_ptr<StronglyConnectedComponents> NewStronglyConnectedComponents(const Graph &G)
    {
        return make_unique<StronglyConnectedComponents>(G);
    }

    inline unique_ptr<Partition> StronglyConnectedComponentsGetPartition(const StronglyConnectedComponents &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline void StronglyConnectedComponentsGetComponentSizes(const StronglyConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void StronglyConnectedComponentsGetComponents(const StronglyConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }

    inline unique_ptr<WeaklyConnectedComponents> NewWeaklyConnectedComponents(const Graph &G)
    {
        return make_unique<WeaklyConnectedComponents>(G);
    }

    inline unique_ptr<Partition> WeaklyConnectedComponentsGetPartition(const WeaklyConnectedComponents &algo)
    {
        return make_unique<Partition>(algo.getPartition());
    }

    inline void WeaklyConnectedComponentsGetComponentSizes(const WeaklyConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        for (auto &&pair : algo.getComponentSizes())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void WeaklyConnectedComponentsGetComponents(const WeaklyConnectedComponents &algo, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        node i = 0;
        for (auto &&block : algo.getComponents())
        {
            for (auto &&n : block)
            {
                ks.push_back(i);
                vs.push_back(n);
            }
            ++i;
        }
    }
}

#endif // NK_COMPONENT_H
