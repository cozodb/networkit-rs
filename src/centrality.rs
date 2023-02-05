use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::{Algorithm, DynAlgorithm},
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
