//
// Created by Ziyang Hu on 2023/1/30.
//

#include "bridge.h"

namespace NetworKit
{
    unique_ptr<vector<edgeweight>> MakeWeightVector(rust::Slice<const edgeweight> wt)
    {
        vector<edgeweight> v{wt.begin(), wt.end()};
        return make_unique<vector<edgeweight>>(v);
    }

    unique_ptr<vector<count>> MakeCountVector(rust::Slice<const count> wt)
    {
        vector<count> v{wt.begin(), wt.end()};
        return make_unique<vector<count>>(v);
    }
}