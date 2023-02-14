use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
    tools::NodeIter,
};


pub trait LinkPredictor {
    fn set_graph(&mut self, g: &crate::Graph);
    fn run(&mut self, u: u64, v: u64);
    fn run_on(&mut self, pairs: &[(u64, u64)]) -> Vec<((u64, u64), f64)>;
    fn run_all(&mut self) -> Vec<((u64, u64), f64)>;
}