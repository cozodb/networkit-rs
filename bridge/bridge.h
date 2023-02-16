//
// Created by Ziyang Hu on 2023/1/29.
//

#ifndef NK_BRIDGE_H
#define NK_BRIDGE_H

#include <graph.h>
#include <community.h>
#include <scd.h>
#include <coarsening.h>
#include <clique.h>
#include <centrality.h>
#include <component.h>
#include <correlation.h>
#include <distance.h>
#include <misc.h>
#include <generators.h>
#include <link_prediction.h>
#include <sparsification.h>
// networkit.dynamics
// networkit.graphio
// networkit.matching
// networkit.randomization
// networkit.simulation

namespace NetworKit
{
    unique_ptr<vector<edgeweight>> MakeWeightVector(rust::Slice<const edgeweight> wt);
    unique_ptr<vector<count>> MakeCountVector(rust::Slice<const count> ct);
}

#endif // NK_BRIDGE_H
