pub mod base;
pub(crate) mod bridge;
pub mod builder;
pub mod centrality;
pub mod clique;
pub mod coarsening;
pub mod community;
pub mod component;
pub mod cover;
pub mod graph;
pub mod partition;
pub mod scd;
pub mod tools;

pub use cover::Cover;
pub use graph::Graph;
pub use partition::Partition;

pub trait QualityMeasure {
    fn get_quality(&mut self, partition: &Partition, graph: &Graph) -> f64;
}
