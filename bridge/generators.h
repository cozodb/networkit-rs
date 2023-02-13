
#ifndef NK_GENERATORS_H
#define NK_GENERATORS_H

#include "rust/cxx.h"
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/generators/ChungLuGenerator.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>

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

    inline unique_ptr<ChungLuGenerator> NewChungLuGenerator(const std::vector<count> &degreeSequence)
    {
        return make_unique<ChungLuGenerator>(degreeSequence);
    }

    inline unique_ptr<Graph> ChungLuGeneratorGenerate(ChungLuGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<ClusteredRandomGraphGenerator> NewClusteredRandomGraphGenerator(count n, count k, double pIntra, double pInter)
    {
        return make_unique<ClusteredRandomGraphGenerator>(n, k, pIntra, pInter);
    }

    inline unique_ptr<Graph> ClusteredRandomGraphGeneratorGenerate(ClusteredRandomGraphGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline unique_ptr<Partition> ClusteredRandomGraphGeneratorGetCommunities(ClusteredRandomGraphGenerator &algo)
    {
        return make_unique<Partition>(algo.getCommunities());
    }

    inline unique_ptr<DorogovtsevMendesGenerator> NewDorogovtsevMendesGenerator(count nNodes)
    {
        return make_unique<DorogovtsevMendesGenerator>(nNodes);
    }

    inline unique_ptr<Graph> DorogovtsevMendesGeneratorGenerate(DorogovtsevMendesGenerator &algo)
    {
        return make_unique<Graph>(algo.generate());
    }

    inline uint8_t convertEventType(GraphEvent::Type t)
    {
        switch (t)
        {
        case GraphEvent::Type::NODE_ADDITION:
            return 0;
            break;
        case GraphEvent::Type::NODE_REMOVAL:
            return 1;
            break;
        case GraphEvent::Type::NODE_RESTORATION:
            return 2;
            break;
        case GraphEvent::Type::EDGE_ADDITION:
            return 3;
            break;
        case GraphEvent::Type::EDGE_REMOVAL:
            return 4;
            break;
        case GraphEvent::Type::EDGE_WEIGHT_UPDATE:
            return 5;
            break;
        case GraphEvent::Type::EDGE_WEIGHT_INCREMENT:
            return 6;
            break;
        case GraphEvent::Type::TIME_STEP:
            return 7;
            break;

        default:
            return 255;
            break;
        }
    }

    inline unique_ptr<DynamicDorogovtsevMendesGenerator> NewDynamicDorogovtsevMendesGenerator()
    {
        return make_unique<DynamicDorogovtsevMendesGenerator>();
    }

    inline void DynamicDorogovtsevMendesGeneratorGenerate(DynamicDorogovtsevMendesGenerator &algo, count nSteps, rust::Vec<uint8_t> &tps, rust::Vec<node> &us, rust::Vec<node> &vs, rust::Vec<edgeweight> &ws)
    {
        for (auto &&evt : algo.generate(nSteps))
        {
            tps.push_back(convertEventType(evt.type));
            us.push_back(evt.u);
            vs.push_back(evt.v);
            ws.push_back(evt.w);
        }
    }
}

#endif // NK_GENERATORS_H