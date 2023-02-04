#ifndef NK_SCD_H
#define NK_SCD_H

#include "rust/cxx.h"
#include <networkit/scd/ApproximatePageRank.hpp>
#include <networkit/scd/CliqueDetect.hpp>
#include <networkit/scd/CombinedSCD.hpp>
#include <networkit/scd/GCE.hpp>
#include <networkit/scd/LFMLocal.hpp>
#include <networkit/scd/LocalT.hpp>
#include <networkit/scd/LocalTightnessExpansion.hpp>
#include <networkit/scd/PageRankNibble.hpp>
#include <networkit/scd/RandomBFS.hpp>
#include <networkit/scd/SCDGroundTruthComparison.hpp>
#include <networkit/scd/SetConductance.hpp>
#include <networkit/scd/TCE.hpp>
#include <networkit/scd/TwoPhaseL.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<ApproximatePageRank> NewApproximatePageRank(const Graph &g, double alpha, double epsilon = 1e-12)
    {
        return make_unique<ApproximatePageRank>(g, alpha, epsilon);
    }

    inline void ApproximatePageRankRun(ApproximatePageRank &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline unique_ptr<CliqueDetect> NewCliqueDetect(const Graph &g)
    {
        return make_unique<CliqueDetect>(g);
    }

    inline void CliqueDetectRun(CliqueDetect &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void CliqueDetectExpandOneCommunity(CliqueDetect &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> CliqueDetectAsBase(unique_ptr<CliqueDetect> algo)
    {
        return algo;
    }

    inline unique_ptr<CombinedSCD> NewCombinedSCD(const Graph &g, SelectiveCommunityDetector &first, SelectiveCommunityDetector &second)
    {
        return make_unique<CombinedSCD>(g, first, second);
    }

    inline void CombinedSCDRun(CombinedSCD &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void CombinedSCDExpandOneCommunity(CombinedSCD &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> CombinedSCDAsBase(unique_ptr<CombinedSCD> algo)
    {
        return algo;
    }

    inline unique_ptr<GCE> NewGCE(const Graph &g, rust::Str Q)
    {
        string q(Q);
        return make_unique<GCE>(g, q);
    }

    inline void GCERun(GCE &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void GCEExpandOneCommunity(GCE &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> GCEAsBase(unique_ptr<GCE> algo)
    {
        return algo;
    }

    inline unique_ptr<LFMLocal> NewLFMLocal(const Graph &g, double alpha = 1.0)
    {
        return make_unique<LFMLocal>(g, alpha);
    }

    inline void LFMLocalRun(LFMLocal &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void LFMLocalExpandOneCommunity(LFMLocal &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> LFMLocalAsBase(unique_ptr<LFMLocal> algo)
    {
        return algo;
    }

    inline unique_ptr<LocalT> NewLocalT(const Graph &g)
    {
        return make_unique<LocalT>(g);
    }

    inline void LocalTRun(LocalT &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void LocalTExpandOneCommunity(LocalT &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> LocalTAsBase(unique_ptr<LocalT> algo)
    {
        return algo;
    }

    inline unique_ptr<LocalTightnessExpansion> NewLocalTightnessExpansion(const Graph &g, double alpha = 1.0)
    {
        return make_unique<LocalTightnessExpansion>(g, alpha);
    }

    inline void LocalTightnessExpansionRun(LocalTightnessExpansion &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void LocalTightnessExpansionExpandOneCommunity(LocalTightnessExpansion &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> LocalTightnessExpansionAsBase(unique_ptr<LocalTightnessExpansion> algo)
    {
        return algo;
    }

    inline unique_ptr<PageRankNibble> NewPageRankNibble(const Graph &g, double alpha, double epsilon)
    {
        return make_unique<PageRankNibble>(g, alpha, epsilon);
    }

    inline void PageRankNibbleRun(PageRankNibble &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void PageRankNibbleExpandOneCommunity(PageRankNibble &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> PageRankNibbleAsBase(unique_ptr<PageRankNibble> algo)
    {
        return algo;
    }

    inline unique_ptr<RandomBFS> NewRandomBFS(const Graph &g, const Cover &c)
    {
        return make_unique<RandomBFS>(g, c);
    }

    inline void RandomBFSRun(RandomBFS &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void RandomBFSExpandOneCommunity(RandomBFS &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> RandomBFSAsBase(unique_ptr<RandomBFS> algo)
    {
        return algo;
    }

    inline unique_ptr<SCDGroundTruthComparison> NewSCDGroundTruthComparison(const Graph &g, const Cover &groundTruth, rust::Slice<const node> ks, rust::Slice<const node> vs, bool ignoreSeeds = false)
    {
        map<node, set<node>> found;
        for (size_t i = 0; i < ks.length(); ++i)
        {
            auto k = ks[i];
            auto v = vs[i];
            if (found.find(k) == found.end())
            {
                set<node> s;
                s.insert(v);
                found[k] = s;
            }
            else
            {
                found[k].insert(v);
            }
        }
        return make_unique<SCDGroundTruthComparison>(g, groundTruth, found, ignoreSeeds);
    }

    inline void SCDGroundTruthComparisonGetIndividualJaccard(const SCDGroundTruthComparison &algo, rust::Vec<index> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.getIndividualJaccard())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void SCDGroundTruthComparisonGetIndividualPrecision(const SCDGroundTruthComparison &algo, rust::Vec<index> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.getIndividualPrecision())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void SCDGroundTruthComparisonGetIndividualRecall(const SCDGroundTruthComparison &algo, rust::Vec<index> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.getIndividualRecall())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline void SCDGroundTruthComparisonGetIndividualF1(const SCDGroundTruthComparison &algo, rust::Vec<index> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.getIndividualF1())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }

    inline unique_ptr<SetConductance> NewSetConductance(const Graph &g, rust::Slice<const node> community)
    {
        set<node> c(community.begin(), community.end());
        return make_unique<SetConductance>(g, c);
    }

    inline unique_ptr<TCE> NewTCE(const Graph &g, bool refine, bool useJaccard)
    {
        return make_unique<TCE>(g, refine, useJaccard);
    }

    inline void TCERun(TCE &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void TCEExpandOneCommunity(TCE &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> TCEAsBase(unique_ptr<TCE> algo)
    {
        return algo;
    }

    inline unique_ptr<TwoPhaseL> NewTwoPhaseL(const Graph &g)
    {
        return make_unique<TwoPhaseL>(g);
    }

    inline void TwoPhaseLRun(TwoPhaseL &algo, rust::Slice<const node> nodes, rust::Vec<node> &ks, rust::Vec<node> &vs)
    {
        set<node> s(nodes.begin(), nodes.end());
        for (auto &&pair : algo.run(s))
        {
            for (auto &&val : pair.second)
            {
                ks.push_back(pair.first);
                vs.push_back(val);
            }
        }
    }
    inline void TwoPhaseLExpandOneCommunity(TwoPhaseL &algo, rust::Slice<const node> seeds, rust::Vec<node> &rs)
    {
        set<node> s(seeds.begin(), seeds.end());
        for (auto &&v : algo.expandOneCommunity(s))
        {
            rs.push_back(v);
        }
    }

    inline unique_ptr<SelectiveCommunityDetector> TwoPhaseLAsBase(unique_ptr<TwoPhaseL> algo)
    {
        return algo;
    }
}

#endif // NK_SCD_H