use cxx::UniquePtr;

use crate::{
    bridge::{self, *},
    tools::NodeIter,
};

pub trait LinkPredictor {
    fn set_graph(&mut self, g: &crate::Graph);
    fn run(&mut self, u: u64, v: u64) -> f64;
    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)>;
    fn run_all(&mut self) -> Vec<((u64, u64), f64)>;
}

pub trait EvaluationMetric {
    fn set_test_graph(&mut self, g: &crate::Graph);
    fn get_curve(
        &mut self,
        prediction: &[((u64, u64), f64)],
        num_threshold: u64,
    ) -> Vec<(f64, f64)>;
    fn get_area_under_curve(&self, curve: &[(f64, f64)]) -> f64;
    fn get_last_area_under_curve(&self) -> f64;
}

pub struct AdamicAdarIndex {
    inner: UniquePtr<bridge::AdamicAdarIndex>,
}

impl AdamicAdarIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewAdamicAdarIndex(g),
        }
    }
}

impl LinkPredictor for AdamicAdarIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        AdamicAdarIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        AdamicAdarIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct AdjustedRandIndex {
    inner: UniquePtr<bridge::AdjustedRandIndex>,
}

impl AdjustedRandIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewAdjustedRandIndex(g),
        }
    }
}

impl LinkPredictor for AdjustedRandIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        AdjustedRandIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        AdjustedRandIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct AlgebraicDistanceIndex {
    inner: UniquePtr<bridge::AlgebraicDistanceIndex>,
}

impl AlgebraicDistanceIndex {
    pub fn new(
        g: &crate::Graph,
        num_systems: u64,
        num_iterations: u64,
        omega: f64,
        norm: u64,
    ) -> Self {
        Self {
            inner: NewAlgebraicDistanceIndex(g, num_systems, num_iterations, omega, norm),
        }
    }
    pub fn preprocess(&mut self) {
        self.inner.pin_mut().preprocess()
    }
}

impl LinkPredictor for AlgebraicDistanceIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        AlgebraicDistanceIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        AlgebraicDistanceIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct CommonNeighborsIndex {
    inner: UniquePtr<bridge::CommonNeighborsIndex>,
}

impl CommonNeighborsIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewCommonNeighborsIndex(g),
        }
    }
}

impl LinkPredictor for CommonNeighborsIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        CommonNeighborsIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        CommonNeighborsIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct ROCMetric {
    inner: UniquePtr<bridge::ROCMetric>,
}

impl ROCMetric {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewROCMetric(g),
        }
    }
}

impl EvaluationMetric for ROCMetric {
    fn set_test_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setTestGraph(g)
    }

    fn get_curve(
        &mut self,
        prediction: &[((u64, u64), f64)],
        num_threshold: u64,
    ) -> Vec<(f64, f64)> {
        let mut us = Vec::with_capacity(prediction.len());
        let mut vs = Vec::with_capacity(prediction.len());
        let mut ws = Vec::with_capacity(prediction.len());
        for ((u, v), w) in prediction {
            us.push(*u);
            vs.push(*v);
            ws.push(*w);
        }
        let mut xs = vec![];
        let mut ys = vec![];
        ROCMetricGetCurve(
            self.inner.pin_mut(),
            &us,
            &vs,
            &ws,
            num_threshold,
            &mut xs,
            &mut ys,
        );
        xs.into_iter().zip(ys).collect()
    }

    fn get_area_under_curve(&self, curve: &[(f64, f64)]) -> f64 {
        let mut xs = Vec::with_capacity(curve.len());
        let mut ys = Vec::with_capacity(curve.len());
        for (x, y) in curve {
            xs.push(*x);
            ys.push(*y);
        }
        ROCMetricGetAreaUnderCurve(&self.inner, &xs, &ys)
    }

    fn get_last_area_under_curve(&self) -> f64 {
        self.inner.getAreaUnderCurve()
    }
}

pub struct PrecisionRecallMetric {
    inner: UniquePtr<bridge::PrecisionRecallMetric>,
}

impl PrecisionRecallMetric {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewPrecisionRecallMetric(g),
        }
    }
}

impl EvaluationMetric for PrecisionRecallMetric {
    fn set_test_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setTestGraph(g)
    }

    fn get_curve(
        &mut self,
        prediction: &[((u64, u64), f64)],
        num_threshold: u64,
    ) -> Vec<(f64, f64)> {
        let mut us = Vec::with_capacity(prediction.len());
        let mut vs = Vec::with_capacity(prediction.len());
        let mut ws = Vec::with_capacity(prediction.len());
        for ((u, v), w) in prediction {
            us.push(*u);
            vs.push(*v);
            ws.push(*w);
        }
        let mut xs = vec![];
        let mut ys = vec![];
        PrecisionRecallMetricGetCurve(
            self.inner.pin_mut(),
            &us,
            &vs,
            &ws,
            num_threshold,
            &mut xs,
            &mut ys,
        );
        xs.into_iter().zip(ys).collect()
    }

    fn get_area_under_curve(&self, curve: &[(f64, f64)]) -> f64 {
        let mut xs = Vec::with_capacity(curve.len());
        let mut ys = Vec::with_capacity(curve.len());
        for (x, y) in curve {
            xs.push(*x);
            ys.push(*y);
        }
        PrecisionRecallMetricGetAreaUnderCurve(&self.inner, &xs, &ys)
    }

    fn get_last_area_under_curve(&self) -> f64 {
        self.inner.getAreaUnderCurve()
    }
}

pub struct JaccardIndex {
    inner: UniquePtr<bridge::JaccardIndex>,
}

impl JaccardIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewJaccardIndex(g),
        }
    }
}

impl LinkPredictor for JaccardIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        JaccardIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        JaccardIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct KatzIndex {
    inner: UniquePtr<bridge::KatzIndex>,
}

impl KatzIndex {
    pub fn new(g: &crate::Graph, max_path_length: Option<u64>, damping_value: Option<f64>) -> Self {
        Self {
            inner: NewKatzIndex(
                g,
                max_path_length.unwrap_or(5),
                damping_value.unwrap_or(0.005),
            ),
        }
    }
}

impl LinkPredictor for KatzIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        KatzIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        KatzIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct LinkThresholder;

impl LinkThresholder {
    pub fn by_count(prediction: &[((u64, u64), f64)], num_links: u64) -> Vec<(u64, u64)> {
        let mut us = Vec::with_capacity(prediction.len());
        let mut vs = Vec::with_capacity(prediction.len());
        let mut ws = Vec::with_capacity(prediction.len());
        for ((u, v), w) in prediction {
            us.push(*u);
            vs.push(*v);
            ws.push(*w);
        }
        let mut src = vec![];
        let mut dst = vec![];
        LinkThresholderByCount(&us, &vs, &ws, num_links, &mut src, &mut dst);
        src.into_iter().zip(dst).collect()
    }
    pub fn by_percentage(
        prediction: &[((u64, u64), f64)],
        percentage_links: f64,
    ) -> Vec<(u64, u64)> {
        let mut us = Vec::with_capacity(prediction.len());
        let mut vs = Vec::with_capacity(prediction.len());
        let mut ws = Vec::with_capacity(prediction.len());
        for ((u, v), w) in prediction {
            us.push(*u);
            vs.push(*v);
            ws.push(*w);
        }
        let mut src = vec![];
        let mut dst = vec![];
        LinkThresholderByPercentage(&us, &vs, &ws, percentage_links, &mut src, &mut dst);
        src.into_iter().zip(dst).collect()
    }
    pub fn by_score(prediction: &[((u64, u64), f64)], min_score: f64) -> Vec<(u64, u64)> {
        let mut us = Vec::with_capacity(prediction.len());
        let mut vs = Vec::with_capacity(prediction.len());
        let mut ws = Vec::with_capacity(prediction.len());
        for ((u, v), w) in prediction {
            us.push(*u);
            vs.push(*v);
            ws.push(*w);
        }
        let mut src = vec![];
        let mut dst = vec![];
        LinkThresholderByScore(&us, &vs, &ws, min_score, &mut src, &mut dst);
        src.into_iter().zip(dst).collect()
    }
}

pub struct MissingLinksFinder {
    inner: UniquePtr<bridge::MissingLinksFinder>,
}

impl MissingLinksFinder {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewMissingLinksFinder(g),
        }
    }
    pub fn find_at_distance(&mut self, k: u64) -> impl Iterator<Item = (u64, u64)> {
        let mut src = vec![];
        let mut dst = vec![];
        MissingLinksFinderFindAtDistance(self.inner.pin_mut(), k, &mut src, &mut dst);
        src.into_iter().zip(dst)
    }
    pub fn find_from_node(&mut self, u: u64, k: u64) -> impl Iterator<Item = (u64, u64)> {
        let mut src = vec![];
        let mut dst = vec![];
        MissingLinksFinderFindFromNode(self.inner.pin_mut(), u, k, &mut src, &mut dst);
        src.into_iter().zip(dst)
    }
}

pub struct NeighborhoodDistanceIndex {
    inner: UniquePtr<bridge::NeighborhoodDistanceIndex>,
}

impl NeighborhoodDistanceIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewNeighborhoodDistanceIndex(g),
        }
    }
}

impl LinkPredictor for NeighborhoodDistanceIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        NeighborhoodDistanceIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        NeighborhoodDistanceIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct NeighbourhoodUtility;

impl NeighbourhoodUtility {
    pub fn get_neighbours_union(g: &crate::Graph, u: u64, v: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: NeighborhoodUtilityGetNeighborsUnion(g, u, v),
        }
    }
    pub fn get_common_neighbours(g: &crate::Graph, u: u64, v: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: NeighborhoodUtilityGetCommonNeighbors(g, u, v),
        }
    }
}

pub struct NeighborsMeasureIndex {
    inner: UniquePtr<bridge::NeighborsMeasureIndex>,
}

impl NeighborsMeasureIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewNeighborsMeasureIndex(g),
        }
    }
}

impl LinkPredictor for NeighborsMeasureIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        NeighborsMeasureIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        NeighborsMeasureIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct PreferentialAttachmentIndex {
    inner: UniquePtr<bridge::PreferentialAttachmentIndex>,
}

impl PreferentialAttachmentIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewPreferentialAttachmentIndex(g),
        }
    }
}

impl LinkPredictor for PreferentialAttachmentIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        PreferentialAttachmentIndexRunOn(
            self.inner.pin_mut(),
            &src,
            &dst,
            &mut ks,
            &mut vs,
            &mut ws,
        );
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        PreferentialAttachmentIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct RandomLinkSampler;

impl RandomLinkSampler {
    pub fn by_count(g: &crate::Graph, num_links: u64) -> crate::Graph {
        RandomLinkSamplerByCount(g, num_links).into()
    }
    pub fn by_percentage(g: &crate::Graph, percentage: f64) -> crate::Graph {
        RandomLinkSamplerByPercentage(g, percentage).into()
    }
}

pub struct ResourceAllocationIndex {
    inner: UniquePtr<bridge::ResourceAllocationIndex>,
}

impl ResourceAllocationIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewResourceAllocationIndex(g),
        }
    }
}

impl LinkPredictor for ResourceAllocationIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        ResourceAllocationIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        ResourceAllocationIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct SameCommunityIndex {
    inner: UniquePtr<bridge::SameCommunityIndex>,
}

impl SameCommunityIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewSameCommunityIndex(g),
        }
    }
}

impl LinkPredictor for SameCommunityIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        SameCommunityIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        SameCommunityIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct TotalNeighborsIndex {
    inner: UniquePtr<bridge::TotalNeighborsIndex>,
}

impl TotalNeighborsIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewTotalNeighborsIndex(g),
        }
    }
}

impl LinkPredictor for TotalNeighborsIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        TotalNeighborsIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        TotalNeighborsIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct UDegreeIndex {
    inner: UniquePtr<bridge::UDegreeIndex>,
}

impl UDegreeIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewUDegreeIndex(g),
        }
    }
}

impl LinkPredictor for UDegreeIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        UDegreeIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        UDegreeIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}

pub struct VDegreeIndex {
    inner: UniquePtr<bridge::VDegreeIndex>,
}

impl VDegreeIndex {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewVDegreeIndex(g),
        }
    }
}

impl LinkPredictor for VDegreeIndex {
    fn set_graph(&mut self, g: &crate::Graph) {
        self.inner.pin_mut().setGraph(g)
    }

    fn run(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().run(u, v)
    }

    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)> {
        let mut src = Vec::with_capacity(pairs.len());
        let mut dst = Vec::with_capacity(pairs.len());
        for (s, d) in pairs {
            src.push(*s);
            dst.push(*d);
        }
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        VDegreeIndexRunOn(self.inner.pin_mut(), &src, &dst, &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }

    fn run_all(&mut self) -> Vec<((u64, u64), f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        VDegreeIndexRunAll(self.inner.pin_mut(), &mut ks, &mut vs, &mut ws);
        ks.into_iter().zip(vs).zip(ws).collect()
    }
}
