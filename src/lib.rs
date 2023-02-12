pub mod base;
pub(crate) mod bridge;
pub mod builder;
pub mod centrality;
pub mod clique;
pub mod coarsening;
pub mod community;
pub mod component;
pub mod correlation;
pub mod cover;
pub mod distance;
pub mod embedding;
pub mod flow;
pub mod graph;
pub mod partition;
pub mod scd;
pub mod tools;
pub mod generators;

pub use cover::Cover;
pub use graph::Graph;
pub use partition::Partition;

pub trait QualityMeasure {
    fn get_quality(&mut self, partition: &Partition, graph: &Graph) -> f64;
}
