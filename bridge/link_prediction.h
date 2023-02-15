
#ifndef NK_LINK_PREDICTION_H
#define NK_LINK_PREDICTION_H

#include "rust/cxx.h"
#include <networkit/linkprediction/AdamicAdarIndex.hpp>
#include <networkit/linkprediction/AdjustedRandIndex.hpp>
#include <networkit/linkprediction/AlgebraicDistanceIndex.hpp>
#include <networkit/linkprediction/CommonNeighborsIndex.hpp>
#include <networkit/linkprediction/JaccardIndex.hpp>
#include <networkit/linkprediction/KatzIndex.hpp>
#include <networkit/linkprediction/ROCMetric.hpp>
#include <networkit/linkprediction/PrecisionRecallMetric.hpp>
#include <networkit/linkprediction/LinkThresholder.hpp>
#include <networkit/linkprediction/MissingLinksFinder.hpp>
#include <networkit/linkprediction/NeighborhoodDistanceIndex.hpp>
#include <networkit/linkprediction/NeighborsMeasureIndex.hpp>
#include <networkit/linkprediction/PreferentialAttachmentIndex.hpp>
#include <networkit/linkprediction/RandomLinkSampler.hpp>
#include <networkit/linkprediction/ResourceAllocationIndex.hpp>
#include <networkit/linkprediction/SameCommunityIndex.hpp>
#include <networkit/linkprediction/TotalNeighborsIndex.hpp>
#include <networkit/linkprediction/UDegreeIndex.hpp>
#include <networkit/linkprediction/VDegreeIndex.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<AdamicAdarIndex> NewAdamicAdarIndex(const Graph &G)
    {
        return make_unique<AdamicAdarIndex>(G);
    }

    inline void AdamicAdarIndexRunOn(AdamicAdarIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void AdamicAdarIndexRunAll(AdamicAdarIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<AdjustedRandIndex> NewAdjustedRandIndex(const Graph &G)
    {
        return make_unique<AdjustedRandIndex>(G);
    }

    inline void AdjustedRandIndexRunOn(AdjustedRandIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void AdjustedRandIndexRunAll(AdjustedRandIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<AlgebraicDistanceIndex> NewAlgebraicDistanceIndex(const Graph &G, count numberSystems, count numberIterations, double omega = 0.5,
                                                                        index norm = 2)
    {
        return make_unique<AlgebraicDistanceIndex>(G, numberSystems, numberIterations, omega, norm);
    }

    inline void AlgebraicDistanceIndexRunOn(AlgebraicDistanceIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void AlgebraicDistanceIndexRunAll(AlgebraicDistanceIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<CommonNeighborsIndex> NewCommonNeighborsIndex(const Graph &G)
    {
        return make_unique<CommonNeighborsIndex>(G);
    }

    inline void CommonNeighborsIndexRunOn(CommonNeighborsIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void CommonNeighborsIndexRunAll(CommonNeighborsIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<PrecisionRecallMetric> NewPrecisionRecallMetric(const Graph &G)
    {
        return make_unique<PrecisionRecallMetric>(G);
    }
    inline void PrecisionRecallMetricGetCurve(PrecisionRecallMetric &algo, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const double> ws, count numThresholds, rust::Vec<double> &xs, rust::Vec<double> &ys)
    {
        vector<LinkPredictor::prediction> predictions{};
        for (size_t i = 0; i < us.length(); ++i)
        {
            predictions.emplace_back(pair(pair(us[i], vs[i]), ws[i]));
        }
        auto res = algo.getCurve(predictions, numThresholds);
        for (size_t i = 0; i < res.first.size(); ++i)
        {
            xs.push_back(res.first[i]);
            ys.push_back(res.second[i]);
        }
    }
    inline double PrecisionRecallMetricGetAreaUnderCurve(const PrecisionRecallMetric &algo, rust::Slice<const double> xs, rust::Slice<const double> ys)
    {
        vector<double> xs_{xs.begin(), xs.end()};
        vector<double> ys_{xs.begin(), ys.end()};

        return algo.getAreaUnderCurve(pair(xs_, ys_));
    }

    inline unique_ptr<ROCMetric> NewROCMetric(const Graph &G)
    {
        return make_unique<ROCMetric>(G);
    }
    inline void ROCMetricGetCurve(ROCMetric &algo, rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const double> ws, count numThresholds, rust::Vec<double> &xs, rust::Vec<double> &ys)
    {
        vector<LinkPredictor::prediction> predictions{};
        for (size_t i = 0; i < us.length(); ++i)
        {
            predictions.emplace_back(pair(pair(us[i], vs[i]), ws[i]));
        }
        auto res = algo.getCurve(predictions, numThresholds);
        for (size_t i = 0; i < res.first.size(); ++i)
        {
            xs.push_back(res.first[i]);
            ys.push_back(res.second[i]);
        }
    }
    inline double ROCMetricGetAreaUnderCurve(const ROCMetric &algo, rust::Slice<const double> xs, rust::Slice<const double> ys)
    {
        vector<double> xs_{xs.begin(), xs.end()};
        vector<double> ys_{xs.begin(), ys.end()};

        return algo.getAreaUnderCurve(pair(xs_, ys_));
    }

    inline unique_ptr<JaccardIndex> NewJaccardIndex(const Graph &G)
    {
        return make_unique<JaccardIndex>(G);
    }

    inline void JaccardIndexRunOn(JaccardIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void JaccardIndexRunAll(JaccardIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<KatzIndex> NewKatzIndex(const Graph &G, count maxPathLength = 5, double dampingValue = 0.005)
    {
        return make_unique<KatzIndex>(G, maxPathLength, dampingValue);
    }

    inline void KatzIndexRunOn(KatzIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void KatzIndexRunAll(KatzIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void LinkThresholderByCount(rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const double> ws, count numLinks, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        vector<LinkPredictor::prediction> predictions{};
        for (size_t i = 0; i < us.length(); ++i)
        {
            predictions.emplace_back(pair(pair(us[i], vs[i]), ws[i]));
        }
        for (auto &&p : LinkThresholder::byCount(predictions, numLinks))
        {
            src.push_back(p.first);
            dst.push_back(p.second);
        }
    }

    inline void LinkThresholderByPercentage(rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const double> ws, double percentageLinks, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        vector<LinkPredictor::prediction> predictions{};
        for (size_t i = 0; i < us.length(); ++i)
        {
            predictions.emplace_back(pair(pair(us[i], vs[i]), ws[i]));
        }
        for (auto &&p : LinkThresholder::byPercentage(predictions, percentageLinks))
        {
            src.push_back(p.first);
            dst.push_back(p.second);
        }
    }

    inline void LinkThresholderByScore(rust::Slice<const node> us, rust::Slice<const node> vs, rust::Slice<const double> ws, double minScore, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        vector<LinkPredictor::prediction> predictions{};
        for (size_t i = 0; i < us.length(); ++i)
        {
            predictions.emplace_back(pair(pair(us[i], vs[i]), ws[i]));
        }
        for (auto &&p : LinkThresholder::byScore(predictions, minScore))
        {
            src.push_back(p.first);
            dst.push_back(p.second);
        }
    }

    inline unique_ptr<MissingLinksFinder> NewMissingLinksFinder(const Graph &G)
    {
        return make_unique<MissingLinksFinder>(G);
    }

    inline void MissingLinksFinderFindAtDistance(MissingLinksFinder &algo, count k, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&p : algo.findAtDistance(k))
        {
            src.push_back(p.first);
            dst.push_back(p.second);
        }
    }
    inline void MissingLinksFinderFindFromNode(MissingLinksFinder &algo, node u, count k, rust::Vec<node> &src, rust::Vec<node> &dst)
    {
        for (auto &&p : algo.findFromNode(u, k))
        {
            src.push_back(p.first);
            dst.push_back(p.second);
        }
    }

    inline unique_ptr<NeighborhoodDistanceIndex> NewNeighborhoodDistanceIndex(const Graph &G)
    {
        return make_unique<NeighborhoodDistanceIndex>(G);
    }

    inline void NeighborhoodDistanceIndexRunOn(NeighborhoodDistanceIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void NeighborhoodDistanceIndexRunAll(NeighborhoodDistanceIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<vector<node>> NeighborhoodUtilityGetNeighborsUnion(const Graph &G, node u, node v)
    {
        return make_unique<vector<node>>(NeighborhoodUtility::getNeighborsUnion(G, u, v));
    }

    inline unique_ptr<vector<node>> NeighborhoodUtilityGetCommonNeighbors(const Graph &G, node u, node v)
    {
        return make_unique<vector<node>>(NeighborhoodUtility::getCommonNeighbors(G, u, v));
    }

    inline unique_ptr<NeighborsMeasureIndex> NewNeighborsMeasureIndex(const Graph &G)
    {
        return make_unique<NeighborsMeasureIndex>(G);
    }

    inline void NeighborsMeasureIndexRunOn(NeighborsMeasureIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void NeighborsMeasureIndexRunAll(NeighborsMeasureIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<PreferentialAttachmentIndex> NewPreferentialAttachmentIndex(const Graph &G)
    {
        return make_unique<PreferentialAttachmentIndex>(G);
    }

    inline void PreferentialAttachmentIndexRunOn(PreferentialAttachmentIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void PreferentialAttachmentIndexRunAll(PreferentialAttachmentIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<Graph> RandomLinkSamplerByCount(const Graph &G, count numLinks)
    {
        return make_unique<Graph>(RandomLinkSampler::byCount(G, numLinks));
    }

    inline unique_ptr<Graph> RandomLinkSamplerByPercentage(const Graph &G, double percentage)
    {
        return make_unique<Graph>(RandomLinkSampler::byPercentage(G, percentage));
    }

    inline unique_ptr<ResourceAllocationIndex> NewResourceAllocationIndex(const Graph &G)
    {
        return make_unique<ResourceAllocationIndex>(G);
    }

    inline void ResourceAllocationIndexRunOn(ResourceAllocationIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void ResourceAllocationIndexRunAll(ResourceAllocationIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<SameCommunityIndex> NewSameCommunityIndex(const Graph &G)
    {
        return make_unique<SameCommunityIndex>(G);
    }

    inline void SameCommunityIndexRunOn(SameCommunityIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void SameCommunityIndexRunAll(SameCommunityIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<TotalNeighborsIndex> NewTotalNeighborsIndex(const Graph &G)
    {
        return make_unique<TotalNeighborsIndex>(G);
    }

    inline void TotalNeighborsIndexRunOn(TotalNeighborsIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void TotalNeighborsIndexRunAll(TotalNeighborsIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<UDegreeIndex> NewUDegreeIndex(const Graph &G)
    {
        return make_unique<UDegreeIndex>(G);
    }

    inline void UDegreeIndexRunOn(UDegreeIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void UDegreeIndexRunAll(UDegreeIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline unique_ptr<VDegreeIndex> NewVDegreeIndex(const Graph &G)
    {
        return make_unique<VDegreeIndex>(G);
    }

    inline void VDegreeIndexRunOn(VDegreeIndex &algo, rust::Slice<const node> src, rust::Slice<const node> dst, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {
        vector<pair<node, node>> pairs{};
        for (size_t i = 0; i < src.length(); ++i)
        {
            pairs.emplace_back(pair(src[i], dst[i]));
        }
        for (auto &&res : algo.runOn(pairs))
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }

    inline void VDegreeIndexRunAll(VDegreeIndex &algo, rust::Vec<node> &ks, rust::Vec<node> &vs, rust::Vec<double> &ws)
    {

        for (auto &&res : algo.runAll())
        {
            ks.push_back(res.first.first);
            vs.push_back(res.first.second);
            ws.push_back(res.second);
        }
    }
}
#endif // NK_LINK_PREDICTION_H