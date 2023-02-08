use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::{Algorithm, DynAlgorithm, EdgeDirection},
    bridge::{self, *},
    community::ValueIter,
    tools::NodeIter,
};

pub trait ConnectedComponents {
    fn extract_largest_connected_component(g: &crate::Graph, compact_graph: bool) -> crate::Graph;
    fn number_of_components(&self) -> u64;
    fn component_of_node(&self, u: u64) -> u64;
    fn get_partition(&self) -> crate::Partition;
    fn get_component_sizes(&self) -> BTreeMap<u64, u64>;
    fn get_components(&self) -> Vec<Vec<u64>>;
}
