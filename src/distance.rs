use std::collections::BTreeMap;

use cxx::{CxxVector, UniquePtr};
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    base::DynAlgorithm,
    bridge::{self, *},
    community::ValueIter,
    tools::NodeIter,
};

pub trait STSP: Algorithm {
    fn get_path(&self) -> NodeIter;
    fn get_predecessors(&self) -> NodeIter;
    fn get_distance(&self) -> f64;
    fn get_distances(&self) -> ValueIter;
    fn set_source(&mut self, u: u64);
    fn set_target(&mut self, v: u64);
    fn set_targets(&mut self, vs: &[u64]);
    fn get_target_index_map(&self) -> BTreeMap<u64, u64>;
}

pub trait SSSP: Algorithm {
    fn distance(&self, t: u64) -> f64;
    fn get_distances(&mut self) -> ValueIter;
    fn number_of_paths(&self, t: u64) -> f64;
    fn get_predecessors(&self, t: u64) -> NodeIter;
    fn get_path(&self, t: u64, forward: bool) -> NodeIter;
    fn get_paths(&self, t: u64, forward: bool) -> Vec<Vec<u64>>;
    fn get_node_sorted_by_distances(&self) -> NodeIter;
    fn get_num_reachable_nodes(&self) -> u64;
    fn set_source(&mut self, u: u64);
    fn set_target(&mut self, v: u64);
    fn get_sum_of_distances(&self) -> f64;
}

pub trait DynSSSP: SSSP + DynAlgorithm {
    fn modified(&mut self) -> bool;
    fn set_target_node(&mut self, t: u64);
}

pub trait NodeDistance {
    fn preprocess(&mut self);
    fn distance(&mut self, u: u64, v: u64) -> f64;
    fn get_edge_scores(&self) -> ValueIter;
}

pub struct APSP {
    inner: UniquePtr<bridge::APSP>,
}

impl APSP {
    pub fn new(g: &crate::Graph) -> Self {
        Self { inner: NewAPSP(g) }
    }
    pub fn get_distance(&self, u: u64, v: u64) -> f64 {
        self.inner.getDistance(u, v)
    }
    pub fn get_distances(&self) -> BTreeMap<u64, BTreeMap<u64, f64>> {
        let mut ret: BTreeMap<u64, BTreeMap<u64, f64>> = BTreeMap::new();
        let mut r = vec![];
        let n = APSPGetDistances(&self.inner, &mut r);
        for (i, wt) in r.into_iter().enumerate() {
            let i = i as u64;
            let src = i / n;
            let dst = i % n;
            ret.entry(src).or_default().insert(dst, wt);
        }
        ret
    }
}

impl Algorithm for APSP {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct AStar {
    inner: UniquePtr<bridge::AStar>,
    _heu: UniquePtr<CxxVector<f64>>,
}

impl AStar {
    pub fn new(g: &crate::Graph, heuristic: &[f64], src: u64, dst: u64, store_pred: bool) -> Self {
        let heu = MakeWeightVector(heuristic);
        let inner = NewAStar(g, &heu, src, dst, store_pred);
        Self { inner, _heu: heu }
    }
}

impl Algorithm for AStar {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl STSP for AStar {
    fn get_path(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: AStarGetPath(&self.inner),
        }
    }

    fn get_predecessors(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: AStarGetPredecessors(&self.inner),
        }
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn set_targets(&mut self, vs: &[u64]) {
        AStarSetTargets(self.inner.pin_mut(), vs);
    }

    fn get_target_index_map(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        AStarGetTargetIndexMap(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }

    fn get_distance(&self) -> f64 {
        self.inner.getDistance()
    }

    fn get_distances(&self) -> ValueIter {
        ValueIter {
            inner: AStarGetDistances(&self.inner),
            at: 0,
        }
    }
}

pub struct AdamicAdarDistance {
    inner: UniquePtr<bridge::AdamicAdarDistance>,
}

impl AdamicAdarDistance {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewAdamicAdarDistance(g),
        }
    }
}

impl NodeDistance for AdamicAdarDistance {
    fn preprocess(&mut self) {
        self.inner.pin_mut().preprocess()
    }

    fn distance(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().distance(u, v)
    }

    fn get_edge_scores(&self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: AdamicAdarDistanceGetEdgeScores(&self.inner),
        }
    }
}

pub struct AlgebraicDistance {
    inner: UniquePtr<bridge::AlgebraicDistance>,
}

impl AlgebraicDistance {
    pub fn new(
        g: &crate::Graph,
        n_systems: Option<u64>,
        n_iterations: Option<u64>,
        omega: Option<f64>,
        norm: Option<u64>,
        with_edge_scores: bool,
    ) -> Self {
        Self {
            inner: NewAlgebraicDistance(
                g,
                n_systems.unwrap_or(10),
                n_iterations.unwrap_or(30),
                omega.unwrap_or(0.5),
                norm.unwrap_or(0),
                with_edge_scores,
            ),
        }
    }
}

impl NodeDistance for AlgebraicDistance {
    fn preprocess(&mut self) {
        self.inner.pin_mut().preprocess()
    }

    fn distance(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().distance(u, v)
    }

    fn get_edge_scores(&self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: AlgebraicDistanceGetEdgeScores(&self.inner),
        }
    }
}

pub struct AllSimplePaths {
    inner: UniquePtr<bridge::AllSimplePaths>,
}

impl AllSimplePaths {
    pub fn new(g: &crate::Graph, src: u64, dst: u64, cutoff: u64) -> Self {
        Self {
            inner: NewAllSimplePaths(g, src, dst, cutoff),
        }
    }

    pub fn number_of_simple_paths(&mut self) -> u64 {
        self.inner.pin_mut().numberOfSimplePaths()
    }

    pub fn get_all_simple_paths(&mut self) -> Vec<Vec<u64>> {
        let mut temp = vec![];
        AllSimplePathsGetAllSimplePaths(self.inner.pin_mut(), &mut temp);
        let mut ret = vec![];
        let mut cur = vec![];
        for n in temp {
            if n == u64::MAX {
                ret.push(cur);
                cur = vec![];
            } else {
                cur.push(n);
            }
        }
        ret
    }
}

impl Algorithm for AllSimplePaths {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct Dijkstra {
    inner: UniquePtr<bridge::Dijkstra>,
}

impl Dijkstra {
    pub fn new(
        g: &crate::Graph,
        src: u64,
        store_paths: bool,
        store_nodes_sorted_by_distance: bool,
        dst: Option<u64>,
    ) -> Self {
        Self {
            inner: NewDijkstra(
                g,
                src,
                store_paths,
                store_nodes_sorted_by_distance,
                dst.unwrap_or(u64::MAX),
            ),
        }
    }
}

impl Algorithm for Dijkstra {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl SSSP for Dijkstra {
    fn distance(&self, t: u64) -> f64 {
        self.inner.distance(t)
    }

    fn get_distances(&mut self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: DijkstraGetDistances(self.inner.pin_mut()),
        }
    }

    fn number_of_paths(&self, t: u64) -> f64 {
        self.inner._numberOfPaths(t)
    }

    fn get_predecessors(&self, t: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: DijkstraGetPredecessors(&self.inner, t),
        }
    }

    fn get_path(&self, t: u64, forward: bool) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: DijkstraGetPath(&self.inner, t, forward),
        }
    }

    fn get_paths(&self, t: u64, forward: bool) -> Vec<Vec<u64>> {
        let mut temp = vec![];
        DijkstraGetPaths(&self.inner, t, forward, &mut temp);
        let mut ret = vec![];
        let mut cur = vec![];
        for n in temp {
            if n == u64::MAX {
                ret.push(cur);
                cur = vec![];
            } else {
                cur.push(n);
            }
        }
        ret
    }

    fn get_node_sorted_by_distances(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: DijkstraGetNodeSortedByDistance(&self.inner),
        }
    }

    fn get_num_reachable_nodes(&self) -> u64 {
        self.inner.getReachableNodes()
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn get_sum_of_distances(&self) -> f64 {
        self.inner.getSumOfDistances()
    }
}

pub struct BidirectionalBFS {
    inner: UniquePtr<bridge::BidirectionalBFS>,
}

impl BidirectionalBFS {
    pub fn new(g: &crate::Graph, src: u64, dst: u64, store_pred: bool) -> Self {
        let inner = NewBidirectionalBFS(g, src, dst, store_pred);
        Self { inner }
    }
}

impl Algorithm for BidirectionalBFS {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl STSP for BidirectionalBFS {
    fn get_path(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BidirectionalBFSGetPath(&self.inner),
        }
    }

    fn get_predecessors(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BidirectionalBFSGetPredecessors(&self.inner),
        }
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn set_targets(&mut self, vs: &[u64]) {
        BidirectionalBFSSetTargets(self.inner.pin_mut(), vs);
    }

    fn get_target_index_map(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        BidirectionalBFSGetTargetIndexMap(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }

    fn get_distance(&self) -> f64 {
        self.inner.getDistance()
    }

    fn get_distances(&self) -> ValueIter {
        ValueIter {
            inner: BidirectionalBFSGetDistances(&self.inner),
            at: 0,
        }
    }
}

pub struct BidirectionalDijkstra {
    inner: UniquePtr<bridge::BidirectionalDijkstra>,
}

impl BidirectionalDijkstra {
    pub fn new(g: &crate::Graph, src: u64, dst: u64, store_pred: bool) -> Self {
        let inner = NewBidirectionalDijkstra(g, src, dst, store_pred);
        Self { inner }
    }
}

impl Algorithm for BidirectionalDijkstra {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl STSP for BidirectionalDijkstra {
    fn get_path(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BidirectionalDijkstraGetPath(&self.inner),
        }
    }

    fn get_predecessors(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BidirectionalDijkstraGetPredecessors(&self.inner),
        }
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn set_targets(&mut self, vs: &[u64]) {
        BidirectionalDijkstraSetTargets(self.inner.pin_mut(), vs);
    }

    fn get_target_index_map(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        BidirectionalDijkstraGetTargetIndexMap(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }

    fn get_distance(&self) -> f64 {
        self.inner.getDistance()
    }

    fn get_distances(&self) -> ValueIter {
        ValueIter {
            inner: BidirectionalDijkstraGetDistances(&self.inner),
            at: 0,
        }
    }
}

pub struct CommuteTimeDistance {
    inner: UniquePtr<bridge::CommuteTimeDistance>,
}

impl CommuteTimeDistance {
    pub fn new(g: &crate::Graph, tol: Option<f64>) -> Self {
        let inner = NewCommuteTimeDistance(g, tol.unwrap_or(0.1));
        Self { inner }
    }

    pub fn run_approximation(&mut self) {
        self.inner.pin_mut().runApproximation()
    }
    pub fn run_parallel_approximation(&mut self) {
        self.inner.pin_mut().runParallelApproximation()
    }
    pub fn get_startup_time(&self) -> u64 {
        self.inner.getSetupTime()
    }
    pub fn run_single_pair(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().runSinglePair(u, v)
    }
    pub fn distance(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().distance(u, v)
    }
    pub fn run_single_source(&mut self, u: u64) -> f64 {
        self.inner.pin_mut().runSingleSource(u)
    }
}

impl Algorithm for CommuteTimeDistance {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct Diameter {
    inner: UniquePtr<bridge::Diameter>,
}

#[derive(Default, Copy, Clone, Eq, PartialEq)]
#[repr(u8)]

pub enum DiameterAlgorithm {
    #[default]
    Automatic = 0,
    Exact = 1,
    EstimatedRange = 2,
    EstimatedSamples = 3,
    EstimatedPedantic = 4,
}

impl Diameter {
    pub fn new(
        g: &crate::Graph,
        algo: DiameterAlgorithm,
        error: Option<f64>,
        n_samples: Option<u64>,
    ) -> Self {
        let inner = NewDiameter(g, algo as u8, error.unwrap_or(-1.), n_samples.unwrap_or(0));
        Self { inner }
    }
    pub fn get_diameter(&self) -> (u64, u64) {
        let mut lower = 0;
        let mut upper = 0;
        DiameterGetDiameter(&self.inner, &mut lower, &mut upper);
        (lower, upper)
    }
}

impl Algorithm for Diameter {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct BFS {
    inner: UniquePtr<bridge::BFS>,
}

impl BFS {
    pub fn new(
        g: &crate::Graph,
        src: u64,
        store_paths: bool,
        store_nodes_sorted_by_distance: bool,
        dst: Option<u64>,
    ) -> Self {
        Self {
            inner: NewBFS(
                g,
                src,
                store_paths,
                store_nodes_sorted_by_distance,
                dst.unwrap_or(u64::MAX),
            ),
        }
    }
}

impl Algorithm for BFS {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl SSSP for BFS {
    fn distance(&self, t: u64) -> f64 {
        self.inner.distance(t)
    }

    fn get_distances(&mut self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: BFSGetDistances(self.inner.pin_mut()),
        }
    }

    fn number_of_paths(&self, t: u64) -> f64 {
        self.inner._numberOfPaths(t)
    }

    fn get_predecessors(&self, t: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BFSGetPredecessors(&self.inner, t),
        }
    }

    fn get_path(&self, t: u64, forward: bool) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BFSGetPath(&self.inner, t, forward),
        }
    }

    fn get_paths(&self, t: u64, forward: bool) -> Vec<Vec<u64>> {
        let mut temp = vec![];
        BFSGetPaths(&self.inner, t, forward, &mut temp);
        let mut ret = vec![];
        let mut cur = vec![];
        for n in temp {
            if n == u64::MAX {
                ret.push(cur);
                cur = vec![];
            } else {
                cur.push(n);
            }
        }
        ret
    }

    fn get_node_sorted_by_distances(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: BFSGetNodeSortedByDistance(&self.inner),
        }
    }

    fn get_num_reachable_nodes(&self) -> u64 {
        self.inner.getReachableNodes()
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn get_sum_of_distances(&self) -> f64 {
        self.inner.getSumOfDistances()
    }
}

pub struct DynAPSP {
    inner: UniquePtr<bridge::DynAPSP>,
}

impl DynAPSP {
    pub fn new(g: &mut crate::Graph) -> Self {
        Self {
            inner: NewDynAPSP(g.inner.pin_mut()),
        }
    }
    pub fn get_distance(&self, u: u64, v: u64) -> f64 {
        self.inner.getDistance(u, v)
    }
    pub fn get_distances(&self) -> BTreeMap<u64, BTreeMap<u64, f64>> {
        let mut ret: BTreeMap<u64, BTreeMap<u64, f64>> = BTreeMap::new();
        let mut r = vec![];
        let n = DynAPSPGetDistances(&self.inner, &mut r);
        for (i, wt) in r.into_iter().enumerate() {
            let i = i as u64;
            let src = i / n;
            let dst = i % n;
            ret.entry(src).or_default().insert(dst, wt);
        }
        ret
    }
}

impl Algorithm for DynAPSP {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl DynAlgorithm for DynAPSP {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynAPSPUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
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
        DynAPSPUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

pub struct DynBFS {
    inner: UniquePtr<bridge::DynBFS>,
}

impl DynBFS {
    pub fn new(g: &crate::Graph, src: u64, store_predecessors: bool) -> Self {
        Self {
            inner: NewDynBFS(g, src, store_predecessors),
        }
    }
}

impl Algorithm for DynBFS {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl SSSP for DynBFS {
    fn distance(&self, t: u64) -> f64 {
        self.inner.distance(t)
    }

    fn get_distances(&mut self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: DynBFSGetDistances(self.inner.pin_mut()),
        }
    }

    fn number_of_paths(&self, t: u64) -> f64 {
        self.inner._numberOfPaths(t)
    }

    fn get_predecessors(&self, t: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: DynBFSGetPredecessors(&self.inner, t),
        }
    }

    fn get_path(&self, t: u64, forward: bool) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: DynBFSGetPath(&self.inner, t, forward),
        }
    }

    fn get_paths(&self, t: u64, forward: bool) -> Vec<Vec<u64>> {
        let mut temp = vec![];
        DynBFSGetPaths(&self.inner, t, forward, &mut temp);
        let mut ret = vec![];
        let mut cur = vec![];
        for n in temp {
            if n == u64::MAX {
                ret.push(cur);
                cur = vec![];
            } else {
                cur.push(n);
            }
        }
        ret
    }

    fn get_node_sorted_by_distances(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: DynBFSGetNodeSortedByDistance(&self.inner),
        }
    }

    fn get_num_reachable_nodes(&self) -> u64 {
        self.inner.getReachableNodes()
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn get_sum_of_distances(&self) -> f64 {
        self.inner.getSumOfDistances()
    }
}

impl DynAlgorithm for DynBFS {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynBFSUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
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
        DynBFSUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}

impl DynSSSP for DynBFS {
    fn modified(&mut self) -> bool {
        self.inner.pin_mut().modified()
    }

    fn set_target_node(&mut self, t: u64) {
        self.inner.pin_mut().setTargetNode(t)
    }
}

pub struct Eccentricity;

impl Eccentricity {
    pub fn get_value(g: &crate::Graph, u: u64) -> (u64, u64) {
        let mut fartherest = 0;
        let mut dist = 0;
        EccentricityGetValue(g, u, &mut fartherest, &mut dist);
        (fartherest, dist)
    }
}

pub struct EffectiveDiameter {
    inner: UniquePtr<bridge::EffectiveDiameter>,
}

impl EffectiveDiameter {
    pub fn new(g: &crate::Graph, ratio: Option<f64>) -> Self {
        let inner = NewEffectiveDiameter(g, ratio.unwrap_or(0.9));
        Self { inner }
    }
    pub fn get_effective_diameter(&self) -> f64 {
        self.inner.getEffectiveDiameter()
    }
}

impl Algorithm for EffectiveDiameter {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct EffectiveDiameterApproximation {
    inner: UniquePtr<bridge::EffectiveDiameterApproximation>,
}

impl EffectiveDiameterApproximation {
    pub fn new(g: &crate::Graph, ratio: Option<f64>, k: Option<u64>, r: Option<u64>) -> Self {
        let inner = NewEffectiveDiameterApproximation(
            g,
            ratio.unwrap_or(0.9),
            k.unwrap_or(64),
            r.unwrap_or(7),
        );
        Self { inner }
    }
    pub fn get_effective_diameter(&self) -> f64 {
        self.inner.getEffectiveDiameter()
    }
}

impl Algorithm for EffectiveDiameterApproximation {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct HopPlotApproximation {
    inner: UniquePtr<bridge::HopPlotApproximation>,
}

impl HopPlotApproximation {
    pub fn new(
        g: &crate::Graph,
        max_distance: Option<u64>,
        k: Option<u64>,
        r: Option<u64>,
    ) -> Self {
        let inner = NewHopPlotApproximation(
            g,
            max_distance.unwrap_or(0),
            k.unwrap_or(64),
            r.unwrap_or(7),
        );
        Self { inner }
    }
    pub fn get_hop_plot(&self) -> BTreeMap<u64, f64> {
        let mut ks = vec![];
        let mut vs = vec![];
        HopPlotApproximationGetHopPlot(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
}

impl Algorithm for HopPlotApproximation {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct JaccardDistance {
    inner: UniquePtr<bridge::JaccardDistance>,
    _triangles: UniquePtr<CxxVector<u64>>,
}

impl JaccardDistance {
    pub fn new(g: &crate::Graph, triangles: &[u64]) -> Self {
        let t = MakeCountVector(triangles);
        let inner = NewJaccardDistance(g, &t);
        Self {
            inner,
            _triangles: t,
        }
    }
}

impl NodeDistance for JaccardDistance {
    fn preprocess(&mut self) {
        self.inner.pin_mut().preprocess()
    }

    fn distance(&mut self, u: u64, v: u64) -> f64 {
        self.inner.pin_mut().distance(u, v)
    }

    fn get_edge_scores(&self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: JaccardDistanceGetEdgeScores(&self.inner),
        }
    }
}

pub struct MultiTargetBFS {
    inner: UniquePtr<bridge::MultiTargetBFS>,
}

impl MultiTargetBFS {
    pub fn new(g: &crate::Graph, src: u64, targets: &[u64]) -> Self {
        let inner = NewMultiTargetBFS(g, src, targets);
        Self { inner }
    }
}

impl Algorithm for MultiTargetBFS {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl STSP for MultiTargetBFS {
    fn get_path(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: MultiTargetBFSGetPath(&self.inner),
        }
    }

    fn get_predecessors(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: MultiTargetBFSGetPredecessors(&self.inner),
        }
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn set_targets(&mut self, vs: &[u64]) {
        MultiTargetBFSSetTargets(self.inner.pin_mut(), vs);
    }

    fn get_target_index_map(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        MultiTargetBFSGetTargetIndexMap(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }

    fn get_distance(&self) -> f64 {
        self.inner.getDistance()
    }

    fn get_distances(&self) -> ValueIter {
        ValueIter {
            inner: MultiTargetBFSGetDistances(&self.inner),
            at: 0,
        }
    }
}

pub struct MultiTargetDijkstra {
    inner: UniquePtr<bridge::MultiTargetDijkstra>,
}

impl MultiTargetDijkstra {
    pub fn new(g: &crate::Graph, src: u64, targets: &[u64]) -> Self {
        let inner = NewMultiTargetDijkstra(g, src, targets);
        Self { inner }
    }
}

impl Algorithm for MultiTargetDijkstra {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl STSP for MultiTargetDijkstra {
    fn get_path(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: MultiTargetDijkstraGetPath(&self.inner),
        }
    }

    fn get_predecessors(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: MultiTargetDijkstraGetPredecessors(&self.inner),
        }
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn set_targets(&mut self, vs: &[u64]) {
        MultiTargetDijkstraSetTargets(self.inner.pin_mut(), vs);
    }

    fn get_target_index_map(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        MultiTargetDijkstraGetTargetIndexMap(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }

    fn get_distance(&self) -> f64 {
        self.inner.getDistance()
    }

    fn get_distances(&self) -> ValueIter {
        ValueIter {
            inner: MultiTargetDijkstraGetDistances(&self.inner),
            at: 0,
        }
    }
}

pub struct NeighborhoodFunction {
    inner: UniquePtr<bridge::NeighborhoodFunction>,
}

impl NeighborhoodFunction {
    pub fn new(g: &crate::Graph) -> Self {
        let inner = NewNeighborhoodFunction(g);
        Self { inner }
    }
    pub fn get_neighbourhood_function(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: NeighborhoodFunctionGetNeighborhoodFunction(&self.inner),
        }
    }
}

impl Algorithm for NeighborhoodFunction {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct NeighborhoodFunctionApproximation {
    inner: UniquePtr<bridge::NeighborhoodFunctionApproximation>,
}

impl NeighborhoodFunctionApproximation {
    pub fn new(g: &crate::Graph, k: Option<u64>, r: Option<u64>) -> Self {
        let inner = NewNeighborhoodFunctionApproximation(g, k.unwrap_or(64), r.unwrap_or(7));
        Self { inner }
    }
    pub fn get_neighbourhood_function(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: NeighborhoodFunctionApproximationGetNeighborhoodFunction(&self.inner),
        }
    }
}

impl Algorithm for NeighborhoodFunctionApproximation {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct NeighborhoodFunctionHeuristic {
    inner: UniquePtr<bridge::NeighborhoodFunctionHeuristic>,
}

#[derive(Default, Copy, Clone, Eq, PartialEq)]
#[repr(u8)]

pub enum NeighborhoodFunctionHeuristicStrategy {
    Random = 0,
    #[default]
    Split = 1,
}

impl NeighborhoodFunctionHeuristic {
    pub fn new(
        g: &crate::Graph,
        n_samples: Option<u64>,
        strategy: NeighborhoodFunctionHeuristicStrategy,
    ) -> Self {
        let inner = NewNeighborhoodFunctionHeuristic(g, n_samples.unwrap_or(0), strategy as u8);
        Self { inner }
    }
    pub fn get_neighbourhood_function(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: NeighborhoodFunctionHeuristicGetNeighborhoodFunction(&self.inner),
        }
    }
}

impl Algorithm for NeighborhoodFunctionHeuristic {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct PrunedLandmarkLabeling {
    inner: UniquePtr<bridge::PrunedLandmarkLabeling>,
}

impl PrunedLandmarkLabeling {
    pub fn new(g: &crate::Graph) -> Self {
        let inner = NewPrunedLandmarkLabeling(g);
        Self { inner }
    }
    pub fn query(&self, u: u64, v: u64) -> u64 {
        self.inner.query(u, v)
    }
}

impl Algorithm for PrunedLandmarkLabeling {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ReverseBFS {
    inner: UniquePtr<bridge::ReverseBFS>,
}

impl ReverseBFS {
    pub fn new(
        g: &crate::Graph,
        src: u64,
        store_paths: bool,
        store_nodes_sorted_by_distance: bool,
        dst: Option<u64>,
    ) -> Self {
        Self {
            inner: NewReverseBFS(
                g,
                src,
                store_paths,
                store_nodes_sorted_by_distance,
                dst.unwrap_or(u64::MAX),
            ),
        }
    }
}

impl Algorithm for ReverseBFS {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl SSSP for ReverseBFS {
    fn distance(&self, t: u64) -> f64 {
        self.inner.distance(t)
    }

    fn get_distances(&mut self) -> ValueIter {
        ValueIter {
            at: 0,
            inner: ReverseBFSGetDistances(self.inner.pin_mut()),
        }
    }

    fn number_of_paths(&self, t: u64) -> f64 {
        self.inner._numberOfPaths(t)
    }

    fn get_predecessors(&self, t: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: ReverseBFSGetPredecessors(&self.inner, t),
        }
    }

    fn get_path(&self, t: u64, forward: bool) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: ReverseBFSGetPath(&self.inner, t, forward),
        }
    }

    fn get_paths(&self, t: u64, forward: bool) -> Vec<Vec<u64>> {
        let mut temp = vec![];
        ReverseBFSGetPaths(&self.inner, t, forward, &mut temp);
        let mut ret = vec![];
        let mut cur = vec![];
        for n in temp {
            if n == u64::MAX {
                ret.push(cur);
                cur = vec![];
            } else {
                cur.push(n);
            }
        }
        ret
    }

    fn get_node_sorted_by_distances(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: ReverseBFSGetNodeSortedByDistance(&self.inner),
        }
    }

    fn get_num_reachable_nodes(&self) -> u64 {
        self.inner.getReachableNodes()
    }

    fn set_source(&mut self, u: u64) {
        self.inner.pin_mut().setSource(u)
    }

    fn set_target(&mut self, v: u64) {
        self.inner.pin_mut().setTarget(v)
    }

    fn get_sum_of_distances(&self) -> f64 {
        self.inner.getSumOfDistances()
    }
}

pub struct Volume;

impl Volume {
    pub fn volume(g: &crate::Graph, radius: f64, n_samples: u64) -> f64 {
        VolumeVolume(g, radius, n_samples)
    }
    pub fn volumes(g: &crate::Graph, radii: &[f64], n_samples: u64) -> ValueIter {
        ValueIter {
            at: 0,
            inner: VolumeVolumes(g, radii, n_samples),
        }
    }
}
