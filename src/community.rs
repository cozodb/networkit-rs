use cxx::{CxxVector, UniquePtr};
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
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
