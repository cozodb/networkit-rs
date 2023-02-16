use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
};

pub trait Sparsifier {
    fn run(&mut self);
    fn get_graph(&mut self) -> crate::Graph;
}