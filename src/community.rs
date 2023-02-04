use cxx::{CxxVector, UniquePtr};
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
    scd::SelectiveCommunityDetectorBase,
    QualityMeasure,
};

pub struct AdjustedRandMeasure {
    inner: UniquePtr<bridge::AdjustedRandMeasure>,
}

impl Default for AdjustedRandMeasure {
    fn default() -> Self {
        Self {
            inner: NewAdjustedRandMeasure(),
        }
    }
}

impl AdjustedRandMeasure {
    pub fn get_dissimilarity(
        &mut self,
        g: &crate::Graph,
        zeta: &crate::Partition,
        eta: &crate::Partition,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarity(g, zeta, eta)
    }
}

pub struct ClusteringGenerator {
    inner: UniquePtr<bridge::ClusteringGenerator>,
}

impl Default for ClusteringGenerator {
    fn default() -> Self {
        Self {
            inner: NewClusteringGenerator(),
        }
    }
}

impl ClusteringGenerator {
    pub fn make_continuous_balanced_clustering(
        &mut self,
        g: &crate::Graph,
        k: u64,
    ) -> crate::Partition {
        CMMakeContinuousBalancedClustering(self.inner.pin_mut(), g, k).into()
    }
    pub fn make_noncontinuous_balanced_clustering(
        &mut self,
        g: &crate::Graph,
        k: u64,
    ) -> crate::Partition {
        CMMakeNoncontinuousBalancedClustering(self.inner.pin_mut(), g, k).into()
    }
    pub fn make_one_clustering(&mut self, g: &crate::Graph) -> crate::Partition {
        CMMakeOneClustering(self.inner.pin_mut(), g).into()
    }
    pub fn make_random_clustering(&mut self, g: &crate::Graph, k: u64) -> crate::Partition {
        CMMakeRandomClustering(self.inner.pin_mut(), g, k).into()
    }
    pub fn make_singleton_clustering(&mut self, g: &crate::Graph) -> crate::Partition {
        CMMakeSingletonClustering(self.inner.pin_mut(), g).into()
    }
}

pub struct ValueIter {
    inner: UniquePtr<CxxVector<f64>>,
    at: usize,
}

impl Iterator for ValueIter {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.get(self.at) {
            Some(v) => {
                self.at += 1;
                Some(*v)
            }
            None => None,
        }
    }
}

pub trait LocalCommunityEvaluation: Algorithm {
    fn get_weighted_average(&self) -> f64;
    fn get_unweighted_average(&self) -> f64;
    fn get_maximum_value(&self) -> f64;
    fn get_minimum_value(&self) -> f64;
    fn get_value(&self, i: u64) -> f64;
    fn get_values(&self) -> ValueIter;
    fn is_small_better(&self) -> bool;
}

pub struct CoverF1Similarity {
    inner: UniquePtr<bridge::CoverF1Similarity>,
}

impl Algorithm for CoverF1Similarity {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for CoverF1Similarity {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: CoverF1SimilarityGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

impl CoverF1Similarity {
    pub fn new(g: &crate::Graph, c: &crate::Cover, reference: &crate::Cover) -> Self {
        Self {
            inner: NewCoverF1Similarity(g, c, reference),
        }
    }
}

pub struct CoverHubDominance {
    inner: UniquePtr<bridge::CoverHubDominance>,
}

impl Algorithm for CoverHubDominance {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for CoverHubDominance {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: CoverHubDominanceGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

impl CoverHubDominance {
    pub fn new(g: &crate::Graph, c: &crate::Cover) -> Self {
        Self {
            inner: NewCoverHubDominance(g, c),
        }
    }
}

pub struct Coverage {
    inner: UniquePtr<bridge::Coverage>,
}

impl Default for Coverage {
    fn default() -> Self {
        Self {
            inner: NewCoverage(),
        }
    }
}

impl QualityMeasure for Coverage {
    fn get_quality(&mut self, partition: &crate::Partition, graph: &crate::Graph) -> f64 {
        self.inner.pin_mut().getQuality(partition, graph)
    }
}

pub trait CommunityDetector: Algorithm {
    fn get_partition(&mut self) -> crate::Partition;
}

pub struct CutClustering {
    inner: UniquePtr<bridge::CutClustering>,
}

impl CutClustering {
    pub fn new(g: &crate::Graph, alpha: f64) -> Self {
        Self {
            inner: NewCutClustering(g, alpha),
        }
    }
    pub fn get_cluster_hierarchy(
        g: &crate::Graph,
    ) -> impl Iterator<Item = (f64, crate::Partition)> {
        let inner = CutClusteringGetClusterHierarchy(g);
        struct It {
            inner: UniquePtr<HierarchyIter>,
        }
        impl Iterator for It {
            type Item = (f64, crate::Partition);

            fn next(&mut self) -> Option<Self::Item> {
                if self.inner.isAtEnd() {
                    None
                } else {
                    let k = self.inner.curKey();
                    let v = self.inner.curVal();
                    self.inner.pin_mut().advance();
                    Some((k, v.into()))
                }
            }
        }
        It { inner }
    }
}

impl CommunityDetector for CutClustering {
    fn get_partition(&mut self) -> crate::Partition {
        CutClusteringGetPartition(self.inner.pin_mut()).into()
    }
}

impl Algorithm for CutClustering {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct EdgeCut {
    inner: UniquePtr<bridge::EdgeCut>,
}

impl Default for EdgeCut {
    fn default() -> Self {
        Self {
            inner: NewEdgeCut(),
        }
    }
}

impl QualityMeasure for EdgeCut {
    fn get_quality(&mut self, partition: &crate::Partition, graph: &crate::Graph) -> f64 {
        self.inner.pin_mut().getQuality(partition, graph)
    }
}

pub mod clustering_tools {
    use crate::bridge::*;

    pub fn communication_graph(g: &crate::Graph, zeta: &mut crate::Partition) -> crate::Graph {
        MakeCommunicationGraph(g, zeta.inner.pin_mut()).into()
    }

    pub fn equal_clusterings(
        zeta: &crate::Partition,
        eta: &crate::Partition,
        g: &mut crate::Graph,
    ) -> bool {
        equalClusterings(zeta, eta, g.inner.pin_mut())
    }

    pub fn get_imbalance(zeta: &crate::Partition) -> f32 {
        getImbalance(zeta)
    }

    pub fn is_one_clustering(g: &crate::Graph, zeta: &crate::Partition) -> bool {
        isOneClustering(g, zeta)
    }
    pub fn is_proper_clustering(g: &crate::Graph, zeta: &crate::Partition) -> bool {
        isProperClustering(g, zeta)
    }
    pub fn is_singleton_clustering(g: &crate::Graph, zeta: &crate::Partition) -> bool {
        isSingletonClustering(g, zeta)
    }
    pub fn weighted_degree_with_cluster(
        g: &crate::Graph,
        zeta: &crate::Partition,
        u: u64,
        cid: u64,
    ) -> u64 {
        weightedDegreeWithCluster(g, zeta, u, cid)
    }
}

pub struct GraphStructuralRandMeasure {
    inner: UniquePtr<bridge::GraphStructuralRandMeasure>,
}

impl Default for GraphStructuralRandMeasure {
    fn default() -> Self {
        Self {
            inner: NewGraphStructuralRandMeasure(),
        }
    }
}

impl GraphStructuralRandMeasure {
    pub fn get_dissimilarity(
        &mut self,
        g: &crate::Graph,
        first: &crate::Partition,
        second: &crate::Partition,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarity(g, first, second)
    }
}

pub struct HubDominance {
    inner: UniquePtr<bridge::HubDominance>,
}

impl Default for HubDominance {
    fn default() -> Self {
        Self {
            inner: NewHubDominance(),
        }
    }
}

impl HubDominance {
    pub fn get_quality(&mut self, partition: &crate::Partition, graph: &crate::Graph) -> f64 {
        self.inner.pin_mut().getQuality(partition, graph)
    }
    pub fn get_quality_for_cover(&mut self, cover: &crate::Cover, graph: &crate::Graph) -> f64 {
        self.inner.pin_mut().getQualityForCover(cover, graph)
    }
}

pub struct IntrapartitionDensity {
    inner: UniquePtr<bridge::IntrapartitionDensity>,
}

impl IntrapartitionDensity {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewIntrapartitionDensity(g, p),
        }
    }
    pub fn get_global(&self) -> f64 {
        self.inner.getGlobal()
    }
}

impl Algorithm for IntrapartitionDensity {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for IntrapartitionDensity {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: IntrapartitionDensityGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

pub struct IsolatedInterpartitionConductance {
    inner: UniquePtr<bridge::IsolatedInterpartitionConductance>,
}

impl IsolatedInterpartitionConductance {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewIsolatedInterpartitionConductance(g, p),
        }
    }
}

impl Algorithm for IsolatedInterpartitionConductance {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for IsolatedInterpartitionConductance {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: IsolatedInterpartitionConductanceGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

pub struct IsolatedInterpartitionExpansion {
    inner: UniquePtr<bridge::IsolatedInterpartitionExpansion>,
}

impl IsolatedInterpartitionExpansion {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewIsolatedInterpartitionExpansion(g, p),
        }
    }
}

impl Algorithm for IsolatedInterpartitionExpansion {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for IsolatedInterpartitionExpansion {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: IsolatedInterpartitionExpansionGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

pub struct JaccardMeasure {
    inner: UniquePtr<bridge::JaccardMeasure>,
}

impl Default for JaccardMeasure {
    fn default() -> Self {
        Self {
            inner: NewJaccardMeasure(),
        }
    }
}

impl JaccardMeasure {
    pub fn get_dissimilarity(
        &mut self,
        g: &crate::Graph,
        zeta: &crate::Partition,
        eta: &crate::Partition,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarity(g, zeta, eta)
    }
}

pub trait OverlappingCommunityDetector: Algorithm {
    fn get_cover(&self) -> crate::Cover;
}

pub struct LFM {
    inner: UniquePtr<bridge::LFM>,
}

impl LFM {
    pub fn new(g: &crate::Graph, scd: &mut SelectiveCommunityDetectorBase) -> Self {
        Self {
            inner: NewLFM(g, scd.inner.pin_mut()),
        }
    }
}

impl OverlappingCommunityDetector for LFM {
    fn get_cover(&self) -> crate::Cover {
        LFMGetCover(&self.inner).into()
    }
}

impl Algorithm for LFM {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct LPDegreeOrdered {
    inner: UniquePtr<bridge::LPDegreeOrdered>,
}

impl LPDegreeOrdered {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewLPDegreeOrdered(g),
        }
    }
    pub fn number_of_iterations(&mut self) -> u64 {
        self.inner.pin_mut().numberOfIterations()
    }
}

impl Algorithm for LPDegreeOrdered {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl CommunityDetector for LPDegreeOrdered {
    fn get_partition(&mut self) -> crate::Partition {
        LPDegreeOrderedGetPartition(self.inner.pin_mut()).into()
    }
}

pub struct LouvainMapEquation {
    inner: UniquePtr<bridge::LouvainMapEquation>,
}

impl LouvainMapEquation {
    pub const RELAX_MAP: &str = "relaxmap";
    pub const SYNCHRONOUS: &str = "synchronous";
    pub const NONE: &str = "none";

    pub fn new(
        g: &crate::Graph,
        hierarchical: bool,
        max_iterations: Option<u64>,
        parallelization_strategy: Option<&str>,
    ) -> Self {
        Self {
            inner: NewLouvainMapEquation(
                g,
                hierarchical,
                max_iterations.unwrap_or(32),
                parallelization_strategy.unwrap_or(Self::RELAX_MAP),
            ),
        }
    }
}

impl Algorithm for LouvainMapEquation {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl CommunityDetector for LouvainMapEquation {
    fn get_partition(&mut self) -> crate::Partition {
        LouvainMapEquationGetPartition(self.inner.pin_mut()).into()
    }
}

pub struct Modularity {
    inner: UniquePtr<bridge::Modularity>,
}

impl Default for Modularity {
    fn default() -> Self {
        Self {
            inner: NewModularity(),
        }
    }
}

impl QualityMeasure for Modularity {
    fn get_quality(&mut self, partition: &crate::Partition, graph: &crate::Graph) -> f64 {
        self.inner.pin_mut().getQuality(partition, graph)
    }
}

pub struct NMIDistance {
    inner: UniquePtr<bridge::NMIDistance>,
}

impl Default for NMIDistance {
    fn default() -> Self {
        Self {
            inner: NewNMIDistance(),
        }
    }
}

impl NMIDistance {
    pub fn get_dissimilarity(
        &mut self,
        g: &crate::Graph,
        zeta: &crate::Partition,
        eta: &crate::Partition,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarity(g, zeta, eta)
    }
}

pub struct NodeStructuralRandMeasure {
    inner: UniquePtr<bridge::NodeStructuralRandMeasure>,
}

impl Default for NodeStructuralRandMeasure {
    fn default() -> Self {
        Self {
            inner: NewNodeStructuralRandMeasure(),
        }
    }
}

impl NodeStructuralRandMeasure {
    pub fn get_dissimilarity(
        &mut self,
        g: &crate::Graph,
        zeta: &crate::Partition,
        eta: &crate::Partition,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarity(g, zeta, eta)
    }
}

pub struct OverlappingNMIDistance {
    inner: UniquePtr<bridge::OverlappingNMIDistance>,
}

#[repr(u8)]
#[derive(Default)]
pub enum OverlappingNMIDistanceNormalization {
    Min = 0,
    GeometricMean = 1,
    ArithmeticMean = 2,
    #[default]
    Max = 3,
    JointEntropy = 4,
}

impl OverlappingNMIDistance {
    pub fn new(normalization: OverlappingNMIDistanceNormalization) -> Self {
        Self {
            inner: NewOverlappingNMIDistance(normalization as u8),
        }
    }

    pub fn get_dissimilarity(
        &mut self,
        g: &crate::Graph,
        zeta: &crate::Partition,
        eta: &crate::Partition,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarity(g, zeta, eta)
    }

    pub fn get_dissimilarity_for_cover(
        &mut self,
        g: &crate::Graph,
        zeta: &crate::Cover,
        eta: &crate::Cover,
    ) -> f64 {
        self.inner.pin_mut().getDissimilarityForCover(g, zeta, eta)
    }
}

pub struct PLM {
    inner: UniquePtr<bridge::PLM>,
}

impl PLM {
    pub const BALANCED: &str = "balanced";
    pub const NONE: &str = "none";
    pub const NONE_RANDOMIZED: &str = "none randomized";
    pub const SIMPLE: &str = "simple";
    pub fn new(
        g: &crate::Graph,
        refine: Option<bool>,
        gamma: Option<f64>,
        par: Option<&str>,
        max_iter: Option<u64>,
        turbo: Option<bool>,
        recurse: Option<bool>,
    ) -> Self {
        Self {
            inner: NewPLM(
                g,
                refine.unwrap_or(false),
                gamma.unwrap_or(1.),
                par.unwrap_or(Self::BALANCED),
                max_iter.unwrap_or(32),
                turbo.unwrap_or(true),
                recurse.unwrap_or(true),
            ),
        }
    }
    pub fn coarsen(g: &crate::Graph, zeta: &crate::Partition) -> (crate::Graph, Vec<u64>) {
        let mut mapping = vec![];
        let ret = PLMCoarsen(g, zeta, &mut mapping);
        (ret.into(), mapping)
    }
    pub fn prolong(
        g: &crate::Graph,
        zeta_coarse: &crate::Partition,
        g_fine: &crate::Graph,
        node_to_meta_node: &[u64],
    ) -> crate::Partition {
        PLMProlong(g, zeta_coarse, g_fine, node_to_meta_node).into()
    }
}

impl Algorithm for PLM {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl CommunityDetector for PLM {
    fn get_partition(&mut self) -> crate::Partition {
        PLMGetPartition(self.inner.pin_mut()).into()
    }
}

pub struct PLP {
    inner: UniquePtr<bridge::PLP>,
}

impl PLP {
    pub fn new(g: &crate::Graph, theta: Option<u64>, max_iterations: Option<u64>) -> Self {
        Self {
            inner: NewPLP(
                g,
                theta.unwrap_or(u64::MAX),
                max_iterations.unwrap_or(u64::MAX),
            ),
        }
    }
    pub fn number_of_iterations(&mut self) -> u64 {
        self.inner.pin_mut().numberOfIterations()
    }
}

impl Algorithm for PLP {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl CommunityDetector for PLP {
    fn get_partition(&mut self) -> crate::Partition {
        PLPGetPartition(self.inner.pin_mut()).into()
    }
}

pub struct ParallelLeiden {
    inner: UniquePtr<bridge::ParallelLeiden>,
}

impl ParallelLeiden {
    pub fn new(
        g: &crate::Graph,
        iterations: Option<u64>,
        randomize: Option<bool>,
        gamma: Option<f64>,
    ) -> Self {
        Self {
            inner: NewParallelLeiden(
                g,
                iterations.unwrap_or(u64::MAX),
                randomize.unwrap_or(true),
                gamma.unwrap_or(1.0),
            ),
        }
    }
}

impl Algorithm for ParallelLeiden {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl CommunityDetector for ParallelLeiden {
    fn get_partition(&mut self) -> crate::Partition {
        ParallelLeidenGetPartition(self.inner.pin_mut()).into()
    }
}

pub struct PartitionFragmentation {
    inner: UniquePtr<bridge::PartitionFragmentation>,
}

impl Algorithm for PartitionFragmentation {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for PartitionFragmentation {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: PartitionFragmentationGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

impl PartitionFragmentation {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewPartitionFragmentation(g, p),
        }
    }
}

pub struct PartitionHubDominance {
    inner: UniquePtr<bridge::PartitionHubDominance>,
}

impl Algorithm for PartitionHubDominance {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for PartitionHubDominance {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: PartitionHubDominanceGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

impl PartitionHubDominance {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewPartitionHubDominance(g, p),
        }
    }
}

pub struct PartitionIntersection {
    inner: UniquePtr<bridge::PartitionIntersection>,
}

impl Default for PartitionIntersection {
    fn default() -> Self {
        Self {
            inner: NewPartitionIntersection(),
        }
    }
}

impl PartitionIntersection {
    pub fn calculate(
        &mut self,
        zeta: &crate::Partition,
        eta: &crate::Partition,
    ) -> crate::Partition {
        PartitionIntersectionCalculate(self.inner.pin_mut(), zeta, eta).into()
    }
}

pub struct StablePartitionNodes {
    inner: UniquePtr<bridge::StablePartitionNodes>,
}

impl Algorithm for StablePartitionNodes {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

impl LocalCommunityEvaluation for StablePartitionNodes {
    fn get_weighted_average(&self) -> f64 {
        self.inner.getWeightedAverage()
    }

    fn get_unweighted_average(&self) -> f64 {
        self.inner.getUnweightedAverage()
    }

    fn get_maximum_value(&self) -> f64 {
        self.inner.getMaximumValue()
    }

    fn get_minimum_value(&self) -> f64 {
        self.inner.getMinimumValue()
    }

    fn get_value(&self, i: u64) -> f64 {
        self.inner.getValue(i)
    }

    fn get_values(&self) -> ValueIter {
        ValueIter {
            inner: StablePartitionNodesGetValues(&self.inner),
            at: 0,
        }
    }

    fn is_small_better(&self) -> bool {
        self.inner.isSmallBetter()
    }
}

impl StablePartitionNodes {
    pub fn new(g: &crate::Graph, p: &crate::Partition) -> Self {
        Self {
            inner: NewStablePartitionNodes(g, p),
        }
    }
    pub fn is_stable(&self, u: u64) -> bool {
        self.inner.isStable(u)
    }
}

// TODO communityGraph: need to implement coarsening module first

// compareCommunities: not implemented in upstream

// TODO: detectCommunities, evalCommunityDetection, inspectCommunities, kCoreCommunityDetection, readCommunities, writeCommunities
