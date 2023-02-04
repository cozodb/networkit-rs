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
    pub fn new(g: &Graph, p: &Partition) -> Self {
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
    pub fn new(g: &Graph, p: &Partition) -> Self {
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
    pub fn new(g: &Graph, p: &Partition) -> Self {
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
    pub fn new(g: &Graph, scd: &mut SelectiveCommunityDetectorBase) -> Self {
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
