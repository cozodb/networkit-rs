use cxx::UniquePtr;

use crate::bridge::{self, *};

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
