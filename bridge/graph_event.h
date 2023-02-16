#ifndef NK_EVENT_H
#define NK_EVENT_H

#include <networkit/base/DynAlgorithm.hpp>

namespace NetworKit
{
    inline GraphEvent toGraphEvent(uint8_t kind, node u, node v, edgeweight ew)
    {
        GraphEvent::Type t;
        switch (kind)
        {
        case 0:
            t = GraphEvent::Type::NODE_ADDITION;
            break;
        case 1:
            t = GraphEvent::Type::NODE_REMOVAL;
            break;
        case 2:
            t = GraphEvent::Type::NODE_RESTORATION;
            break;
        case 3:
            t = GraphEvent::Type::EDGE_ADDITION;
            break;
        case 4:
            t = GraphEvent::Type::EDGE_WEIGHT_UPDATE;
            break;
        case 5:
            t = GraphEvent::Type::EDGE_WEIGHT_INCREMENT;
            break;
        case 6:
            t = GraphEvent::Type::TIME_STEP;
            break;
        default:
            t = GraphEvent::Type::TIME_STEP;
            break;
        }

        return GraphEvent(t, u, v, ew);
    }
}

#endif // NK_EVENT_H