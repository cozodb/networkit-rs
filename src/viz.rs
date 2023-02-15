use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
};

trait GraphLayoutAlgorithm {
    fn run(&mut self);
    fn get_coordinates(&self) -> Vec<(f64, f64)>;
    fn num_edge_crossings(&self) -> u64;
    
}