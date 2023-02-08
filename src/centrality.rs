use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::{Algorithm, DynAlgorithm, EdgeDirection},
    bridge::{self, *},
    community::ValueIter,
    tools::NodeIter,
};

pub trait Centrality: Algorithm {
    fn centralization(&mut self) -> f64;
    fn maximum(&mut self) -> f64;
    fn ranking(&mut self) -> RankIter;
    fn score(&mut self, v: u64) -> f64;
    fn scores(&mut self) -> ValueIter;
}

pub struct RankIter {
    pub(crate) ks: Vec<u64>,
    pub(crate) vs: Vec<f64>,
    pub(crate) at: usize,
}

impl Iterator for RankIter {
    type Item = (u64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.at >= self.ks.len() {
            None
        } else {
            let k = self.ks[self.at];
            let v = self.vs[self.at];
            self.at += 1;
            Some((k, v))
        }
    }
}

pub struct ApproxBetweenness {
    inner: UniquePtr<bridge::ApproxBetweenness>,
}

impl ApproxBetweenness {
    pub fn new(
        g: &crate::Graph,
        epsilon: Option<f64>,
        delta: Option<f64>,
        universal_constant: Option<f64>,
    ) -> Self {
        Self {
            inner: NewApproxBetweenness(
                g,
                epsilon.unwrap_or(0.01),
                delta.unwrap_or(0.1),
                universal_constant.unwrap_or(1.),
            ),
        }
    }
}

impl Centrality for ApproxBetweenness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproxBetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxBetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ApproxBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ApproxCloseness {
    inner: UniquePtr<bridge::ApproxCloseness>,
}

#[derive(Default)]
#[repr(u8)]
pub enum ApproxClosenessType {
    #[default]
    OutBound = 0,
    InBound = 1,
    Sum = 2,
}

impl ApproxCloseness {
    pub fn new(
        g: &crate::Graph,
        n_samples: u64,
        epsilon: Option<f64>,
        normalized: bool,
        t: ApproxClosenessType,
    ) -> Self {
        Self {
            inner: NewApproxCloseness(g, n_samples, epsilon.unwrap_or(0.1), normalized, t as u8),
        }
    }

    pub fn get_square_error_estimates(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxClosenessGetSquareErrorEstimates(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Centrality for ApproxCloseness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproxClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxClosenessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ApproxCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ApproxElectricalCloseness {
    inner: UniquePtr<bridge::ApproxElectricalCloseness>,
}

impl ApproxElectricalCloseness {
    pub fn new(g: &crate::Graph, epsilon: Option<f64>, kappa: Option<f64>) -> Self {
        Self {
            inner: NewApproxElectricalCloseness(g, epsilon.unwrap_or(0.1), kappa.unwrap_or(0.3)),
        }
    }

    pub fn compute_exact_diagonal(&self, tol: Option<f64>) -> ValueIter {
        ValueIter {
            inner: ApproxElectricalClosenessComputeExactDiagonal(&self.inner, tol.unwrap_or(1e-9)),
            at: 0,
        }
    }
    pub fn get_diagonal(&self) -> ValueIter {
        ValueIter {
            inner: ApproxElectricalClosenessGetDiagonal(&self.inner),
            at: 0,
        }
    }
}

impl Centrality for ApproxElectricalCloseness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproxElectricalClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxElectricalClosenessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ApproxElectricalCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ApproxGroupBetweenness {
    inner: UniquePtr<bridge::ApproxGroupBetweenness>,
}

impl ApproxGroupBetweenness {
    pub fn new(g: &crate::Graph, group_size: u64, epsilon: f64) -> Self {
        Self {
            inner: NewApproxGroupBetweenness(g, group_size, epsilon),
        }
    }
    pub fn group_max_betweenness(&self) -> impl Iterator<Item = u64> {
        NodeIter {
            nodes: ApproxGroupBetweennessGroupMaxBetweenness(&self.inner),
            at: 0,
        }
    }
    pub fn score_of_group(&self, nodes: &[u64], normalized: bool) -> impl Iterator<Item = u64> {
        NodeIter {
            nodes: ApproxGroupBetweennessScoreOfGroup(&self.inner, nodes, normalized),
            at: 0,
        }
    }
}

impl Algorithm for ApproxGroupBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ApproxSpanningEdge {
    inner: UniquePtr<bridge::ApproxSpanningEdge>,
}

impl ApproxSpanningEdge {
    pub fn new(g: &crate::Graph, epsilon: Option<f64>) -> Self {
        Self {
            inner: NewApproxSpanningEdge(g, epsilon.unwrap_or(0.1)),
        }
    }

    pub fn scores(&self) -> ValueIter {
        ValueIter {
            inner: ApproxSpanningEdgeScores(&self.inner),
            at: 0,
        }
    }
}

impl Algorithm for ApproxSpanningEdge {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct Betweenness {
    inner: UniquePtr<bridge::Betweenness>,
}

impl Betweenness {
    pub fn new(g: &crate::Graph, normalized: bool, compute_edge_centrality: bool) -> Self {
        Self {
            inner: NewBetweenness(g, normalized, compute_edge_centrality),
        }
    }

    pub fn edge_scores(&mut self) -> ValueIter {
        ValueIter {
            inner: BetweennessEdgeScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Centrality for Betweenness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        BetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: BetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for Betweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct Closeness {
    inner: UniquePtr<bridge::Closeness>,
}

#[derive(Default)]
#[repr(u8)]
pub enum ClosenessVariant {
    #[default]
    Standard = 0,
    Generalized = 1,
}

impl Closeness {
    pub fn new(g: &crate::Graph, normalized: bool, variant: ClosenessVariant) -> Self {
        Self {
            inner: NewCloseness(g, normalized, variant as u8),
        }
    }
}

impl Centrality for Closeness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ClosenessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for Closeness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct CoreDecomposition {
    inner: UniquePtr<bridge::CoreDecomposition>,
}

impl CoreDecomposition {
    pub fn new(
        g: &crate::Graph,
        normalized: bool,
        enforce_bucket_queue_algorithm: bool,
        store_node_order: bool,
    ) -> Self {
        Self {
            inner: NewCoreDecomposition(
                g,
                normalized,
                enforce_bucket_queue_algorithm,
                store_node_order,
            ),
        }
    }

    pub fn get_cover(&self) -> crate::Cover {
        CoreDecompositionGetCover(&self.inner).into()
    }
    pub fn get_partition(&self) -> crate::Partition {
        CoreDecompositionGetPartition(&self.inner).into()
    }
    pub fn get_node_order(&self) -> impl Iterator<Item = u64> {
        NodeIter {
            nodes: CoreDecompositionGetNodeOrder(&self.inner),
            at: 0,
        }
    }
    pub fn max_core_number(&self) -> u64 {
        self.inner.maxCoreNumber()
    }
}

impl Centrality for CoreDecomposition {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        CoreDecompositionRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: CoreDecompositionScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for CoreDecomposition {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct DegreeCentrality {
    inner: UniquePtr<bridge::DegreeCentrality>,
}

impl DegreeCentrality {
    pub fn new(g: &crate::Graph, normalized: bool, out_deg: bool, ignore_self_loops: bool) -> Self {
        Self {
            inner: NewDegreeCentrality(g, normalized, out_deg, ignore_self_loops),
        }
    }
}

impl Centrality for DegreeCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        DegreeCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: DegreeCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for DegreeCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct DynApproxBetweenness {
    inner: UniquePtr<bridge::DynApproxBetweenness>,
}

impl DynApproxBetweenness {
    pub fn new(
        g: &crate::Graph,
        epsilon: Option<f64>,
        delta: Option<f64>,
        store_predecessors: bool,
        universal_constant: Option<f64>,
    ) -> Self {
        Self {
            inner: NewDynApproxBetweenness(
                g,
                epsilon.unwrap_or(0.01),
                delta.unwrap_or(0.1),
                store_predecessors,
                universal_constant.unwrap_or(1.),
            ),
        }
    }
    pub fn get_number_of_samples(&self) -> u64 {
        self.inner.getNumberOfSamples()
    }
}

impl Centrality for DynApproxBetweenness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        DynApproxBetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: DynApproxBetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for DynApproxBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl DynAlgorithm for DynApproxBetweenness {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynApproxBetweennessUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynApproxBetweennessUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

pub struct DynBetweenness {
    inner: UniquePtr<bridge::DynBetweenness>,
}

impl DynBetweenness {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewDynBetweenness(g),
        }
    }
}

impl Centrality for DynBetweenness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        DynBetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: DynBetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for DynBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl DynAlgorithm for DynBetweenness {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynBetweennessUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynBetweennessUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

pub struct DynBetweennessOneNode {
    inner: UniquePtr<bridge::DynBetweennessOneNode>,
}

impl DynBetweennessOneNode {
    pub fn new(g: &mut crate::Graph, x: u64) -> Self {
        Self {
            inner: NewDynBetweennessOneNode(g.inner.pin_mut(), x),
        }
    }
    pub fn run(&mut self) {
        self.inner.pin_mut().run()
    }
    pub fn compute_score(&mut self, e: crate::base::GraphEvent) -> f64 {
        DynBetweennessOneNodeComputeScore(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }
    pub fn get_distance(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().getDistance(u, v)
    }
    pub fn get_sigma(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().getSigma(u, v)
    }
    pub fn get_sigma_x(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().getSigmax(u, v)
    }
    pub fn get_bcx(&mut self) -> f64 {
        self.inner.pin_mut().getbcx()
    }
}

impl DynAlgorithm for DynBetweennessOneNode {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynBetweennessOneNodeUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynBetweennessOneNodeUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

pub struct DynKatzCentrality {
    inner: UniquePtr<bridge::DynKatzCentrality>,
}

impl DynKatzCentrality {
    pub fn new(g: &crate::Graph, k: u64, group_only: bool, tolerance: Option<f64>) -> Self {
        Self {
            inner: NewDynKatzCentrality(g, k, group_only, tolerance.unwrap_or(1e-9)),
        }
    }
    pub fn are_distinguished(&mut self, u: u64, v: u64) -> bool {
        self.inner.pin_mut().areDistinguished(u, v)
    }

    pub fn bound(&mut self, u: u64) -> f64 {
        self.inner.pin_mut().bound(u)
    }

    pub fn top(&mut self, n: Option<u64>) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: DynKatzCentralityTop(self.inner.pin_mut(), n.unwrap_or(0)),
        }
    }
}

impl Centrality for DynKatzCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        DynKatzCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: DynKatzCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for DynKatzCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl DynAlgorithm for DynKatzCentrality {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynKatzCentralityUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynKatzCentralityUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

pub struct DynTopHarmonicCloseness {
    inner: UniquePtr<bridge::DynTopHarmonicCloseness>,
}

impl DynTopHarmonicCloseness {
    pub fn new(g: &crate::Graph, k: u64, use_bfs_bound: bool) -> Self {
        Self {
            inner: NewDynTopHarmonicCloseness(g, k, use_bfs_bound),
        }
    }

    pub fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        DynTopHarmonicClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    pub fn reset(&mut self) {
        self.inner.pin_mut().reset()
    }

    pub fn top_k_nodes_list(&mut self, include_trail: bool) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: DynTopHarmonicClosenessTopkNodesList(self.inner.pin_mut(), include_trail),
        }
    }
    pub fn top_k_scores_list(&mut self, include_trail: bool) -> impl Iterator<Item = f64> {
        ValueIter {
            at: 0,
            inner: DynTopHarmonicClosenessTopkScoresList(self.inner.pin_mut(), include_trail),
        }
    }
}

impl Algorithm for DynTopHarmonicCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl DynAlgorithm for DynTopHarmonicCloseness {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynTopHarmonicClosenessUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynTopHarmonicClosenessUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

pub struct EigenvectorCentrality {
    inner: UniquePtr<bridge::EigenvectorCentrality>,
}

impl EigenvectorCentrality {
    pub fn new(g: &crate::Graph, tol: Option<f64>) -> Self {
        Self {
            inner: NewEigenvectorCentrality(g, tol.unwrap_or(1e-9)),
        }
    }
}

impl Centrality for EigenvectorCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        EigenvectorCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: EigenvectorCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for EigenvectorCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct EstimateBetweenness {
    inner: UniquePtr<bridge::EstimateBetweenness>,
}

impl EstimateBetweenness {
    pub fn new(g: &crate::Graph, n_samples: u64, normalized: bool, parallel: bool) -> Self {
        Self {
            inner: NewEstimateBetweenness(g, n_samples, normalized, parallel),
        }
    }
}

impl Centrality for EstimateBetweenness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        EstimateBetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: EstimateBetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for EstimateBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ForestCentrality {
    inner: UniquePtr<bridge::ForestCentrality>,
}

impl ForestCentrality {
    pub fn new(g: &crate::Graph, root: u64, epsilon: Option<f64>, kappa: Option<f64>) -> Self {
        Self {
            inner: NewForestCentrality(g, root, epsilon.unwrap_or(0.1), kappa.unwrap_or(0.3)),
        }
    }

    pub fn get_number_of_samples(&self) -> u64 {
        self.inner.getNumberOfSamples()
    }

    pub fn get_diagonal(&self) -> impl Iterator<Item = f64> {
        ValueIter {
            inner: ForestCentralityGetDiagonal(&self.inner),
            at: 0,
        }
    }
}

impl Centrality for ForestCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ForestCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ForestCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ForestCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GedWalk {
    inner: UniquePtr<bridge::GedWalk>,
}

#[derive(Default, Copy, Clone)]
#[repr(u8)]
pub enum GedWalkBoundStrategy {
    No,
    Spectral,
    #[default]
    Geometric,
    AdaptiveGeometric,
}

#[derive(Default, Copy, Clone)]
#[repr(u8)]
pub enum GedWalkGreedyStrategy {
    #[default]
    Lazy,
    Stochastic,
}

impl GedWalk {
    pub fn new(
        g: &crate::Graph,
        k: Option<u64>,
        init_epsilon: Option<f64>,
        alpha: Option<f64>,
        bs: GedWalkBoundStrategy,
        gs: GedWalkGreedyStrategy,
        spectral_delta: Option<f64>,
    ) -> Self {
        Self {
            inner: NewGedWalk(
                g,
                k.unwrap_or(1),
                init_epsilon.unwrap_or(0.1),
                alpha.unwrap_or(1.),
                bs as u8,
                gs as u8,
                spectral_delta.unwrap_or(0.5),
            ),
        }
    }
    pub fn get_approximate_score(&self) -> f64 {
        self.inner.getApproximateScore()
    }
    pub fn group_max_ged_walk(&self) -> impl Iterator<Item = u64> {
        NodeIter {
            nodes: GedWalkGroupMaxGedWalk(&self.inner),
            at: 0,
        }
    }
    pub fn score_of_group(&mut self, group: &[u64], epsilon: Option<f64>) -> f64 {
        GedWalkScoreOfGroup(self.inner.pin_mut(), group, epsilon.unwrap_or(0.1))
    }
}

impl Algorithm for GedWalk {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GroupCloseness {
    inner: UniquePtr<bridge::GroupCloseness>,
}

impl GroupCloseness {
    pub fn new(g: &crate::Graph, k: Option<u64>, h: Option<u64>) -> Self {
        Self {
            inner: NewGroupCloseness(g, k.unwrap_or(1), h.unwrap_or(0)),
        }
    }
    pub fn score_of_group(&mut self, group: &[u64]) -> f64 {
        GroupClosenessScoreOfGroup(self.inner.pin_mut(), group)
    }
    pub fn group_max_closeness(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: GroupClosenessGroupMaxCloseness(self.inner.pin_mut()),
        }
    }
    pub fn compute_farness(&self, group: &[u64], h: Option<u64>) -> f64 {
        GroupClosenessComputeFarness(&self.inner, group, h.unwrap_or(u64::MAX))
    }
}

impl Algorithm for GroupCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GroupClosenessGrowShrink {
    inner: UniquePtr<bridge::GroupClosenessGrowShrink>,
}

impl GroupClosenessGrowShrink {
    pub fn new(
        g: &crate::Graph,
        group: &[u64],
        extended: bool,
        insertions: Option<u64>,
        max_iterations: Option<u64>,
    ) -> Self {
        Self {
            inner: NewGroupClosenessGrowShrink(
                g,
                group,
                extended,
                insertions.unwrap_or(0),
                max_iterations.unwrap_or(100),
            ),
        }
    }
    pub fn number_of_iterations(&self) -> u64 {
        self.inner.numberOfIterations()
    }
    pub fn group_max_closeness(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: GroupClosenessGrowShrinkGroupMaxCloseness(&self.inner),
        }
    }
}

impl Algorithm for GroupClosenessGrowShrink {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GroupClosenessLocalSearch {
    inner: UniquePtr<bridge::GroupClosenessLocalSearch>,
}

impl GroupClosenessLocalSearch {
    pub fn new(
        g: &crate::Graph,
        group: &[u64],
        run_grow_shrink: bool,
        max_iterations: Option<u64>,
    ) -> Self {
        Self {
            inner: NewGroupClosenessLocalSearch(
                g,
                group,
                run_grow_shrink,
                max_iterations.unwrap_or(u64::MAX),
            ),
        }
    }
    pub fn number_of_iterations(&self) -> u64 {
        self.inner.numberOfIterations()
    }
    pub fn group_max_closeness(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: GroupClosenessLocalSearchGroupMaxCloseness(&self.inner),
        }
    }
}

impl Algorithm for GroupClosenessLocalSearch {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GroupClosenessLocalSwaps {
    inner: UniquePtr<bridge::GroupClosenessLocalSwaps>,
}

impl GroupClosenessLocalSwaps {
    pub fn new(g: &crate::Graph, group: &[u64], max_swaps: Option<u64>) -> Self {
        Self {
            inner: NewGroupClosenessLocalSwaps(g, group, max_swaps.unwrap_or(100)),
        }
    }
    pub fn number_of_swaps(&self) -> u64 {
        self.inner.numberOfSwaps()
    }
    pub fn group_max_closeness(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: GroupClosenessLocalSwapsGroupMaxCloseness(&self.inner),
        }
    }
}

impl Algorithm for GroupClosenessLocalSwaps {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GroupDegree {
    inner: UniquePtr<bridge::GroupDegree>,
}

impl GroupDegree {
    pub fn new(g: &crate::Graph, k: Option<u64>, count_group_nodes: bool) -> Self {
        Self {
            inner: NewGroupDegree(g, k.unwrap_or(1), count_group_nodes),
        }
    }
    pub fn get_score(&mut self) -> u64 {
        self.inner.pin_mut().getScore()
    }
    pub fn score_of_group(&self, group: &[u64]) -> f64 {
        GroupDegreeScoreOfGroup(&self.inner, group)
    }
    pub fn group_max_degree(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: GroupDegreeGroupMaxDegree(self.inner.pin_mut()),
        }
    }
}

impl Algorithm for GroupDegree {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct GroupHarmonicCloseness {
    inner: UniquePtr<bridge::GroupHarmonicCloseness>,
}

impl GroupHarmonicCloseness {
    pub fn new(g: &crate::Graph, k: Option<u64>) -> Self {
        Self {
            inner: NewGroupHarmonicCloseness(g, k.unwrap_or(1)),
        }
    }
    pub fn score_of_group(g: &crate::Graph, group: &[u64]) -> f64 {
        GroupHarmonicClosenessScoreOfGroup(g, group)
    }
    pub fn group_max_harmonic_degree(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: GroupHarmonicClosenessGroupMaxHarmonicCloseness(self.inner.pin_mut()),
        }
    }
}

impl Algorithm for GroupHarmonicCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct HarmonicCloseness {
    inner: UniquePtr<bridge::HarmonicCloseness>,
}

impl HarmonicCloseness {
    pub fn new(g: &crate::Graph, normalized: bool) -> Self {
        Self {
            inner: NewHarmonicCloseness(g, normalized),
        }
    }
}

impl Centrality for HarmonicCloseness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        HarmonicClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: HarmonicClosenessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for HarmonicCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct KPathCentrality {
    inner: UniquePtr<bridge::KPathCentrality>,
}

impl KPathCentrality {
    pub fn new(g: &crate::Graph, alpha: Option<f64>, k: Option<u64>) -> Self {
        Self {
            inner: NewKPathCentrality(g, alpha.unwrap_or(0.2), k.unwrap_or(0)),
        }
    }
}

impl Centrality for KPathCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        KPathCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: KPathCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for KPathCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct KadabraBetweenness {
    inner: UniquePtr<bridge::KadabraBetweenness>,
}

impl KadabraBetweenness {
    pub fn new(
        g: &crate::Graph,
        err: Option<f64>,
        delta: Option<f64>,
        deterministic: bool,
        k: Option<u64>,
        union_sample: Option<u64>,
        start_factor: Option<u64>,
    ) -> Self {
        Self {
            inner: NewKadabraBetweenness(
                g,
                err.unwrap_or(0.01),
                delta.unwrap_or(0.1),
                deterministic,
                k.unwrap_or(0),
                union_sample.unwrap_or(0),
                start_factor.unwrap_or(100),
            ),
        }
    }

    pub fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        KadabraBetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    pub fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: KadabraBetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }

    pub fn get_number_of_iterations(&self) -> u64 {
        self.inner.getNumberOfIterations()
    }

    pub fn get_omega(&self) -> f64 {
        self.inner.getOmega()
    }

    pub fn top_k_nodes_list(&mut self) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: KadabraBetweennessTopkNodesList(self.inner.pin_mut()),
        }
    }
    pub fn top_k_scores_list(&mut self) -> impl Iterator<Item = f64> {
        ValueIter {
            at: 0,
            inner: KadabraBetweennessTopkScoresList(self.inner.pin_mut()),
        }
    }
}

impl Algorithm for KadabraBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct KatzCentrality {
    inner: UniquePtr<bridge::KatzCentrality>,
}

impl KatzCentrality {
    pub fn new(g: &crate::Graph, alpha: Option<f64>, beta: Option<f64>, tol: Option<f64>) -> Self {
        Self {
            inner: NewKatzCentrality(
                g,
                alpha.unwrap_or(0.),
                beta.unwrap_or(0.2),
                tol.unwrap_or(1e-8),
            ),
        }
    }
    pub fn set_edge_direction(&mut self, dir: EdgeDirection) {
        KatzCentralitySetEdgeDirection(self.inner.pin_mut(), dir == EdgeDirection::OutEdges)
    }
}

impl Centrality for KatzCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        KatzCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: KatzCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for KatzCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct LaplacianCentrality {
    inner: UniquePtr<bridge::LaplacianCentrality>,
}

impl LaplacianCentrality {
    pub fn new(g: &crate::Graph, normalized: bool) -> Self {
        Self {
            inner: NewLaplacianCentrality(g, normalized),
        }
    }
}

impl Centrality for LaplacianCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        LaplacianCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: LaplacianCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for LaplacianCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct LocalClusteringCoefficient {
    inner: UniquePtr<bridge::LocalClusteringCoefficient>,
}

impl LocalClusteringCoefficient {
    pub fn new(g: &crate::Graph, turbo: bool) -> Self {
        Self {
            inner: NewLocalClusteringCoefficient(g, turbo),
        }
    }
}

impl Centrality for LocalClusteringCoefficient {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        LocalClusteringCoefficientRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: LocalClusteringCoefficientScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for LocalClusteringCoefficient {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct LocalPartitionCoverage {
    inner: UniquePtr<bridge::LocalPartitionCoverage>,
}

impl LocalPartitionCoverage {
    pub fn new(g: &crate::Graph, partition: &crate::Partition) -> Self {
        Self {
            inner: NewLocalPartitionCoverage(g, partition),
        }
    }
}

impl Centrality for LocalPartitionCoverage {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        LocalPartitionCoverageRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: LocalPartitionCoverageScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for LocalPartitionCoverage {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct LocalSquareClusteringCoefficient {
    inner: UniquePtr<bridge::LocalSquareClusteringCoefficient>,
}

impl LocalSquareClusteringCoefficient {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewLocalSquareClusteringCoefficient(g),
        }
    }
}

impl Centrality for LocalSquareClusteringCoefficient {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        LocalSquareClusteringCoefficientRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: LocalSquareClusteringCoefficientScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for LocalSquareClusteringCoefficient {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct PageRank {
    inner: UniquePtr<bridge::PageRank>,
}

#[derive(Default)]
#[repr(u8)]
pub enum PageRankNorm {
    L1 = 0,
    #[default]
    L2 = 1,
}

impl PageRank {
    pub fn new(
        g: &crate::Graph,
        damp: Option<f64>,
        tol: Option<f64>,
        normalized: bool,
        distribute_sinks: bool,
    ) -> Self {
        Self {
            inner: NewPageRank(
                g,
                damp.unwrap_or(0.85),
                tol.unwrap_or(1e-9),
                normalized,
                distribute_sinks,
            ),
        }
    }
    pub fn set_max_iterations(&mut self, max_iter: u64) {
        PageRankSetMaxIterations(self.inner.pin_mut(), max_iter)
    }
    pub fn set_norm(&mut self, norm: PageRankNorm) {
        PageRankSetNorm(self.inner.pin_mut(), norm as u8)
    }
    pub fn number_of_iterations(&self) -> u64 {
        self.inner.numberOfIterations()
    }
}

impl Centrality for PageRank {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        PageRankRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: PageRankScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for PageRank {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct PermanenceCentrality {
    inner: UniquePtr<bridge::PermanenceCentrality>,
}

impl PermanenceCentrality {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewPermanenceCentrality(g, p),
        }
    }
    pub fn get_permanence(&mut self, u: u64) -> f64 {
        self.inner.pin_mut().getPermanence(u)
    }
    pub fn get_intra_clustering(&mut self, u: u64) -> f64 {
        self.inner.pin_mut().getIntraClustering(u)
    }
}

impl Algorithm for PermanenceCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct Sfigality {
    inner: UniquePtr<bridge::Sfigality>,
}

impl Sfigality {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewSfigality(g),
        }
    }
}

impl Centrality for Sfigality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        SfigalityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: SfigalityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for Sfigality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct SpanningEdgeCentrality {
    inner: UniquePtr<bridge::SpanningEdgeCentrality>,
}

impl SpanningEdgeCentrality {
    pub fn new(g: &crate::Graph, tol: Option<f64>) -> Self {
        Self {
            inner: NewSpanningEdgeCentrality(g, tol.unwrap_or(0.1)),
        }
    }
    pub fn run_approximation(&mut self) {
        self.inner.pin_mut().runApproximation()
    }
    pub fn run_parallel_approximation(&mut self) {
        self.inner.pin_mut().runParallelApproximation()
    }
}

impl Centrality for SpanningEdgeCentrality {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        SpanningEdgeCentralityRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: SpanningEdgeCentralityScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for SpanningEdgeCentrality {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct TopCloseness {
    inner: UniquePtr<bridge::TopCloseness>,
}

impl TopCloseness {
    pub fn new(g: &crate::Graph, k: u64, first_heu: bool, sec_heu: bool) -> Self {
        Self {
            inner: NewTopCloseness(g, k, first_heu, sec_heu),
        }
    }

    pub fn top_k_nodes_list(&mut self, include_trail: bool) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: TopClosenessTopkNodesList(self.inner.pin_mut(), include_trail),
        }
    }
    pub fn top_k_scores_list(&mut self, include_trail: bool) -> impl Iterator<Item = f64> {
        ValueIter {
            at: 0,
            inner: TopClosenessTopkScoresList(self.inner.pin_mut(), include_trail),
        }
    }
    pub fn restrict_top_k_computation_to_nodes(&mut self, nodes: &[u64]) {
        TopClosenessRestrictTopKComputationToNodes(self.inner.pin_mut(), nodes)
    }
}

impl Algorithm for TopCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct TopHarmonicCloseness {
    inner: UniquePtr<bridge::TopHarmonicCloseness>,
}

impl TopHarmonicCloseness {
    pub fn new(g: &crate::Graph, k: u64, use_nb_bound: bool) -> Self {
        Self {
            inner: NewTopHarmonicCloseness(g, k, use_nb_bound),
        }
    }

    pub fn top_k_nodes_list(&mut self, include_trail: bool) -> impl Iterator<Item = u64> {
        NodeIter {
            at: 0,
            nodes: TopHarmonicClosenessTopkNodesList(self.inner.pin_mut(), include_trail),
        }
    }
    pub fn top_k_scores_list(&mut self, include_trail: bool) -> impl Iterator<Item = f64> {
        ValueIter {
            at: 0,
            inner: TopHarmonicClosenessTopkScoresList(self.inner.pin_mut(), include_trail),
        }
    }
    pub fn restrict_top_k_computation_to_nodes(&mut self, nodes: &[u64]) {
        TopHarmonicClosenessRestrictTopKComputationToNodes(self.inner.pin_mut(), nodes)
    }
}

impl Algorithm for TopHarmonicCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}
