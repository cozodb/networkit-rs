use cxx::{CxxVector, UniquePtr};
use miette::IntoDiagnostic;

use crate::{
    base::{Algorithm, GraphEvent, GraphEventType},
    bridge::{self, *},
    tools::NodeIter,
};

pub trait StaticGraphGenerator {
    fn generate(&mut self) -> crate::Graph;
}

pub trait StaticDegreeSequenceGenerator: StaticGraphGenerator {
    fn is_realizable(&mut self) -> bool;
    fn get_realizable(&self) -> bool;
}

pub trait DynamicGraphGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent>;
}

pub struct BarabasiAlbertGenerator {
    inner: UniquePtr<bridge::BarabasiAlbertGenerator>,
}

impl BarabasiAlbertGenerator {
    pub fn new(k: u64, n_max: u64, n0: Option<u64>, batagelj: bool) -> Self {
        Self {
            inner: NewBarabasiAlbertGenerator(k, n_max, n0.unwrap_or(0), batagelj),
        }
    }
}

impl StaticGraphGenerator for BarabasiAlbertGenerator {
    fn generate(&mut self) -> crate::Graph {
        BarabasiAlbertGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct ChungLuGenerator {
    inner: UniquePtr<bridge::ChungLuGenerator>,
    _seq: UniquePtr<CxxVector<u64>>,
}

impl ChungLuGenerator {
    pub fn new(degree_sequence: &[u64]) -> Self {
        let seq = MakeCountVector(degree_sequence);
        let inner = NewChungLuGenerator(&seq);
        Self { inner, _seq: seq }
    }
}

impl StaticGraphGenerator for ChungLuGenerator {
    fn generate(&mut self) -> crate::Graph {
        ChungLuGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct ClusteredRandomGraphGenerator {
    inner: UniquePtr<bridge::ClusteredRandomGraphGenerator>,
}

impl ClusteredRandomGraphGenerator {
    pub fn new(n: u64, k: u64, p_intra: f64, p_inter: f64) -> Self {
        let inner = NewClusteredRandomGraphGenerator(n, k, p_intra, p_inter);
        Self { inner }
    }
    pub fn get_communities(&mut self) -> crate::Partition {
        ClusteredRandomGraphGeneratorGetCommunities(self.inner.pin_mut()).into()
    }
}

impl StaticGraphGenerator for ClusteredRandomGraphGenerator {
    fn generate(&mut self) -> crate::Graph {
        ClusteredRandomGraphGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct DorogovtsevMendesGenerator {
    inner: UniquePtr<bridge::DorogovtsevMendesGenerator>,
}

impl DorogovtsevMendesGenerator {
    pub fn new(n_nodes: u64) -> Self {
        let inner = NewDorogovtsevMendesGenerator(n_nodes);
        Self { inner }
    }
}

impl StaticGraphGenerator for DorogovtsevMendesGenerator {
    fn generate(&mut self) -> crate::Graph {
        DorogovtsevMendesGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct DynamicDorogovtsevMendesGenerator {
    inner: UniquePtr<bridge::DynamicDorogovtsevMendesGenerator>,
}

impl DynamicDorogovtsevMendesGenerator {
    pub fn new() -> Self {
        let inner = NewDynamicDorogovtsevMendesGenerator();
        Self { inner }
    }
}

impl DynamicGraphGenerator for DynamicDorogovtsevMendesGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent> {
        let mut typs = vec![];
        let mut us = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        DynamicDorogovtsevMendesGeneratorGenerate(
            self.inner.pin_mut(),
            n_steps,
            &mut typs,
            &mut us,
            &mut vs,
            &mut ws,
        );
        typs.into_iter()
            .zip(us)
            .zip(vs)
            .zip(ws)
            .map(|(((t, u), v), w)| GraphEvent {
                kind: GraphEventType::from(t),
                u,
                v,
                ew: w,
            })
            .collect()
    }
}

pub struct DynamicForestFireGenerator {
    inner: UniquePtr<bridge::DynamicForestFireGenerator>,
}

impl DynamicForestFireGenerator {
    pub fn new(p: f64, directed: bool, r: Option<f64>) -> Self {
        let inner = NewDynamicForestFireGenerator(p, directed, r.unwrap_or(1.));
        Self { inner }
    }
}

impl DynamicGraphGenerator for DynamicForestFireGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent> {
        let mut typs = vec![];
        let mut us = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        DynamicForestFireGeneratorGenerate(
            self.inner.pin_mut(),
            n_steps,
            &mut typs,
            &mut us,
            &mut vs,
            &mut ws,
        );
        typs.into_iter()
            .zip(us)
            .zip(vs)
            .zip(ws)
            .map(|(((t, u), v), w)| GraphEvent {
                kind: GraphEventType::from(t),
                u,
                v,
                ew: w,
            })
            .collect()
    }
}

pub struct DynamicHyperbolicGenerator {
    inner: UniquePtr<bridge::DynamicHyperbolicGenerator>,
}

impl DynamicHyperbolicGenerator {
    pub fn new(
        n: u64,
        avg_degree: f64,
        exp: f64,
        t: f64,
        move_each_step: f64,
        move_distance: f64,
    ) -> Self {
        let inner =
            NewDynamicHyperbolicGenerator(n, avg_degree, exp, t, move_each_step, move_distance);
        Self { inner }
    }
    pub fn get_graph(&self) -> crate::Graph {
        DynamicHyperbolicGeneratorGetGraph(&self.inner).into()
    }
    pub fn get_coordinates(&self) -> impl Iterator<Item = (f64, f64)> {
        let mut xs = vec![];
        let mut ys = vec![];
        DynamicHyperbolicGeneratorGetCoordinates(&self.inner, &mut xs, &mut ys);
        xs.into_iter().zip(ys)
    }
}

impl DynamicGraphGenerator for DynamicHyperbolicGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent> {
        let mut typs = vec![];
        let mut us = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        DynamicHyperbolicGeneratorGenerate(
            self.inner.pin_mut(),
            n_steps,
            &mut typs,
            &mut us,
            &mut vs,
            &mut ws,
        );
        typs.into_iter()
            .zip(us)
            .zip(vs)
            .zip(ws)
            .map(|(((t, u), v), w)| GraphEvent {
                kind: GraphEventType::from(t),
                u,
                v,
                ew: w,
            })
            .collect()
    }
}

pub struct DynamicPathGenerator {
    inner: UniquePtr<bridge::DynamicPathGenerator>,
}

impl DynamicPathGenerator {
    pub fn new() -> Self {
        let inner = NewDynamicPathGenerator();
        Self { inner }
    }
}

impl DynamicGraphGenerator for DynamicPathGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent> {
        let mut typs = vec![];
        let mut us = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        DynamicPathGeneratorGenerate(
            self.inner.pin_mut(),
            n_steps,
            &mut typs,
            &mut us,
            &mut vs,
            &mut ws,
        );
        typs.into_iter()
            .zip(us)
            .zip(vs)
            .zip(ws)
            .map(|(((t, u), v), w)| GraphEvent {
                kind: GraphEventType::from(t),
                u,
                v,
                ew: w,
            })
            .collect()
    }
}

pub struct DynamicPubWebGenerator {
    inner: UniquePtr<bridge::DynamicPubWebGenerator>,
}

impl DynamicPubWebGenerator {
    pub fn new(
        num_nodes: u64,
        num_dense_areas: u64,
        neighbourhood_radius: f64,
        max_num_neighbours: u64,
        write_initial_graph_to_stream: bool,
    ) -> Self {
        let inner = NewDynamicPubWebGenerator(
            num_nodes,
            num_dense_areas,
            neighbourhood_radius,
            max_num_neighbours,
            write_initial_graph_to_stream,
        );
        Self { inner }
    }
    pub fn get_graph(&self) -> crate::Graph {
        DynamicPubWebGeneratorGetGraph(&self.inner).into()
    }
    pub fn get_coordinates(&self) -> impl Iterator<Item = (f64, f64)> {
        let mut xs = vec![];
        let mut ys = vec![];
        DynamicPubWebGeneratorGetCoordinates(&self.inner, &mut xs, &mut ys);
        xs.into_iter().zip(ys)
    }
    pub fn get_new_coordinates(&self) -> impl Iterator<Item = (u64, (f64, f64))> {
        let mut ns = vec![];
        let mut xs = vec![];
        let mut ys = vec![];
        DynamicPubWebGeneratorGetNewCoordinates(&self.inner, &mut ns, &mut xs, &mut ys);
        ns.into_iter().zip(xs.into_iter().zip(ys))
    }
}

impl DynamicGraphGenerator for DynamicPubWebGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent> {
        let mut typs = vec![];
        let mut us = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        DynamicPubWebGeneratorGenerate(
            self.inner.pin_mut(),
            n_steps,
            &mut typs,
            &mut us,
            &mut vs,
            &mut ws,
        );
        typs.into_iter()
            .zip(us)
            .zip(vs)
            .zip(ws)
            .map(|(((t, u), v), w)| GraphEvent {
                kind: GraphEventType::from(t),
                u,
                v,
                ew: w,
            })
            .collect()
    }
}

pub struct EdgeSwitchingMarkovChainGenerator {
    inner: UniquePtr<bridge::EdgeSwitchingMarkovChainGenerator>,
    _seq: UniquePtr<CxxVector<u64>>,
}

impl EdgeSwitchingMarkovChainGenerator {
    pub fn new(
        sequence: &[u64],
        ignore_if_not_realizable: bool,
        num_switches_per_edge: u64,
    ) -> Self {
        let seq = MakeCountVector(sequence);
        let inner = NewEdgeSwitchingMarkovChainGenerator(
            &seq,
            ignore_if_not_realizable,
            num_switches_per_edge,
        );
        Self { inner, _seq: seq }
    }
}

impl StaticGraphGenerator for EdgeSwitchingMarkovChainGenerator {
    fn generate(&mut self) -> crate::Graph {
        EdgeSwitchingMarkovChainGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

impl StaticDegreeSequenceGenerator for EdgeSwitchingMarkovChainGenerator {
    fn is_realizable(&mut self) -> bool {
        self.inner.pin_mut().isRealizable()
    }

    fn get_realizable(&self) -> bool {
        self.inner.getRealizable()
    }
}

pub struct ErdosRenyiGenerator {
    inner: UniquePtr<bridge::ErdosRenyiGenerator>,
}

impl ErdosRenyiGenerator {
    pub fn new(n_nodes: u64, prob: f64, directed: bool, self_loops: bool) -> Self {
        let inner = NewErdosRenyiGenerator(n_nodes, prob, directed, self_loops);
        Self { inner }
    }
}

impl StaticGraphGenerator for ErdosRenyiGenerator {
    fn generate(&mut self) -> crate::Graph {
        ErdosRenyiGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct HavelHakimiGenerator {
    inner: UniquePtr<bridge::HavelHakimiGenerator>,
    _seq: UniquePtr<CxxVector<u64>>,
}

impl HavelHakimiGenerator {
    pub fn new(sequence: &[u64], ignore_if_not_realizable: bool) -> Self {
        let seq = MakeCountVector(sequence);
        let inner = NewHavelHakimiGenerator(&seq, ignore_if_not_realizable);
        Self { inner, _seq: seq }
    }
}

impl StaticGraphGenerator for HavelHakimiGenerator {
    fn generate(&mut self) -> crate::Graph {
        HavelHakimiGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

impl StaticDegreeSequenceGenerator for HavelHakimiGenerator {
    fn is_realizable(&mut self) -> bool {
        self.inner.pin_mut().isRealizable()
    }

    fn get_realizable(&self) -> bool {
        self.inner.getRealizable()
    }
}

pub struct HyperbolicGenerator {
    inner: UniquePtr<bridge::HyperbolicGenerator>,
}

impl HyperbolicGenerator {
    pub fn new(n: u64, avg_degree: f64, exp: f64, t: f64) -> Self {
        let inner = NewHyperbolicGenerator(n, avg_degree, exp, t);
        Self { inner }
    }
    pub fn generate_advanced(
        &mut self,
        angles: &[f64],
        radii: &[f64],
        r: f64,
        t: f64,
    ) -> crate::Graph {
        HyperbolicGeneratorGenerateAdvanced(self.inner.pin_mut(), angles, radii, r, t).into()
    }

    pub fn set_leaf_capacity(&mut self, capacity: u64) {
        self.inner.pin_mut().setLeafCapacity(capacity)
    }
    pub fn set_theoretical_split(&mut self, split: bool) {
        self.inner.pin_mut().setTheoreticalSplit(split)
    }
    pub fn set_balance(&mut self, balance: f64) {
        self.inner.pin_mut().setBalance(balance)
    }
}

impl StaticGraphGenerator for HyperbolicGenerator {
    fn generate(&mut self) -> crate::Graph {
        HyperbolicGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct LFRGenerator {
    inner: UniquePtr<bridge::LFRGenerator>,
    _seq: UniquePtr<CxxVector<f64>>,
}

impl LFRGenerator {
    pub fn new(n: u64) -> Self {
        let seq = MakeWeightVector(&[]);
        let inner = NewLFRGenerator(n);
        Self { inner, _seq: seq }
    }
    pub fn set_mu_array(&mut self, mu: &[f64]) {
        let seq = MakeWeightVector(mu);
        self.inner.pin_mut().setMuArray(&seq);
        self._seq = seq;
    }
    pub fn set_mut(&mut self, mu: f64) {
        self.inner.pin_mut().setMu(mu)
    }
    pub fn set_mut_with_binomial_distribution(&mut self, mu: f64) {
        self.inner.pin_mut().setMuWithBinomialDistribution(mu)
    }
    pub fn set_degree_sequence(&mut self, seq: &[u64]) {
        LFRGeneratorSetDegreeSequence(self.inner.pin_mut(), seq)
    }
    pub fn set_community_size_sequence(&mut self, seq: &[u64]) {
        LFRGeneratorSetCommunitySizeSequence(self.inner.pin_mut(), seq)
    }
    pub fn set_partition(&mut self, p: crate::Partition) {
        LFRGeneratorSetPartition(self.inner.pin_mut(), p.inner)
    }
    pub fn generate_powerlaw_degree_sequence(
        &mut self,
        avg_degree: u64,
        max_degree: u64,
        node_degree_exp: f64,
    ) {
        self.inner
            .pin_mut()
            .generatePowerlawDegreeSequence(avg_degree, max_degree, node_degree_exp)
    }
    pub fn generate_powerlaw_community_size_sequence(
        &mut self,
        min_community_size: u64,
        max_community_size: u64,
        community_size_exp: f64,
    ) {
        self.inner.pin_mut().generatePowerlawCommunitySizeSequence(
            min_community_size,
            max_community_size,
            community_size_exp,
        )
    }
    pub fn get_graph(&self) -> crate::Graph {
        LFRGeneratorGetGraph(&self.inner).into()
    }
    pub fn get_partition(&self) -> crate::Partition {
        LFRGeneratorGetPartition(&self.inner).into()
    }
}

impl StaticGraphGenerator for LFRGenerator {
    fn generate(&mut self) -> crate::Graph {
        LFRGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

impl Algorithm for LFRGenerator {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct MocnikGenerator {
    inner: UniquePtr<bridge::MocnikGenerator>,
}

impl MocnikGenerator {
    pub fn new(dim: u64, n: u64, k: f64, weighted: bool) -> Self {
        let inner = NewMocnikGenerator(dim, n, k, weighted);
        Self { inner }
    }
}

impl StaticGraphGenerator for MocnikGenerator {
    fn generate(&mut self) -> crate::Graph {
        MocnikGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct MocnikGeneratorBasic {
    inner: UniquePtr<bridge::MocnikGeneratorBasic>,
}

impl MocnikGeneratorBasic {
    pub fn new(dim: u64, n: u64, k: f64) -> Self {
        let inner = NewMocnikGeneratorBasic(dim, n, k);
        Self { inner }
    }
}

impl StaticGraphGenerator for MocnikGeneratorBasic {
    fn generate(&mut self) -> crate::Graph {
        MocnikGeneratorBasicGenerate(self.inner.pin_mut()).into()
    }
}

pub struct PowerlawDegreeSequence {
    inner: UniquePtr<bridge::PowerlawDegreeSequence>,
}

impl PowerlawDegreeSequence {
    pub fn new(min_deg: u64, max_deg: u64, gamma: f64) -> Self {
        let inner = NewPowerlawDegreeSequence(min_deg, max_deg, gamma);
        Self { inner }
    }
    pub fn get_degree_sequence(&self, num_nodes: u64) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: PowerlawDegreeSequenceGetDegreeSequence(&self.inner, num_nodes),
        }
    }
    pub fn set_minimum_from_average_degree(&mut self, avg_deg: f64) {
        self.inner.pin_mut().setMinimumFromAverageDegree(avg_deg)
    }
    pub fn set_gamma_from_average_degree(&mut self, avg_deg: f64, min_gamma: f64, max_gamma: f64) {
        self.inner
            .pin_mut()
            .setGammaFromAverageDegree(avg_deg, min_gamma, max_gamma)
    }
    pub fn set_minimum_degree(&mut self, min_deg: u64) {
        self.inner.pin_mut().setMinimumDegree(min_deg)
    }
    pub fn get_minimum_degree(&mut self) -> u64 {
        self.inner.getMinimumDegree()
    }
    pub fn get_maximum_degree(&mut self) -> u64 {
        self.inner.getMaximumDegree()
    }
    pub fn set_gamma(&mut self, gamma: f64) {
        self.inner.pin_mut().setGamma(gamma)
    }
    pub fn get_gamma(&self) -> f64 {
        self.inner.getGamma()
    }
    pub fn get_degree(&self) -> u64 {
        self.inner.getDegree()
    }
    pub fn get_expected_average_degree(&self) -> f64 {
        self.inner.getExpectedAverageDegree()
    }
}

impl Algorithm for PowerlawDegreeSequence {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct PubWebGenerator {
    inner: UniquePtr<bridge::PubWebGenerator>,
}

impl PubWebGenerator {
    pub fn new(
        num_nodes: u64,
        num_dense_areas: u64,
        neighbourhood_radius: f64,
        max_num_neighbours: u64,
    ) -> Self {
        let inner = NewPubWebGenerator(
            num_nodes,
            num_dense_areas,
            neighbourhood_radius,
            max_num_neighbours,
        );
        Self { inner }
    }
    pub fn get_coordinates(&self) -> impl Iterator<Item = (f64, f64)> {
        let mut xs = vec![];
        let mut ys = vec![];
        PubWebGeneratorGetCoordinates(&self.inner, &mut xs, &mut ys);
        xs.into_iter().zip(ys)
    }
}

impl StaticGraphGenerator for PubWebGenerator {
    fn generate(&mut self) -> crate::Graph {
        PubWebGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct RegularRingLatticeGenerator {
    inner: UniquePtr<bridge::RegularRingLatticeGenerator>,
}

impl RegularRingLatticeGenerator {
    pub fn new(num_nodes: u64, num_neighbours: u64) -> Self {
        let inner = NewRegularRingLatticeGenerator(num_nodes, num_neighbours);
        Self { inner }
    }
}

impl StaticGraphGenerator for RegularRingLatticeGenerator {
    fn generate(&mut self) -> crate::Graph {
        RegularRingLatticeGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct RmatGenerator {
    inner: UniquePtr<bridge::RmatGenerator>,
}

impl RmatGenerator {
    pub fn new(
        scale: u64,
        edge_factor: u64,
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        weighted: bool,
        reduce_nodes: u64,
    ) -> Self {
        let inner = NewRmatGenerator(scale, edge_factor, a, b, c, d, weighted, reduce_nodes);
        Self { inner }
    }
}

impl StaticGraphGenerator for RmatGenerator {
    fn generate(&mut self) -> crate::Graph {
        RmatGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct WattsStrogatzGenerator {
    inner: UniquePtr<bridge::WattsStrogatzGenerator>,
}

impl WattsStrogatzGenerator {
    pub fn new(n_nodes: u64, n_neighbours: u64, p: f64) -> Self {
        let inner = NewWattsStrogatzGenerator(n_nodes, n_neighbours, p);
        Self { inner }
    }
}

impl StaticGraphGenerator for WattsStrogatzGenerator {
    fn generate(&mut self) -> crate::Graph {
        WattsStrogatzGeneratorGenerate(self.inner.pin_mut()).into()
    }
}
