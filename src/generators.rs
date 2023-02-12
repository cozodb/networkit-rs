use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
    tools::NodeIter,
};

pub trait StaticGraphGenerator {
    fn generate(&mut self) -> crate::Graph;
}


pub struct BarabasiAlbertGenerator {
    inner: UniquePtr<bridge::BarabasiAlbertGenerator>,
}

impl BarabasiAlbertGenerator {
    pub fn new(k: u64, n_max: u64, n0: Option<u64>, batagelj: bool) -> Self {
        Self {
            inner: NewBarabasiAlbertGenerator(k, n_max, n0.unwrap_or(0), batagelj),
        }
    }
}

impl StaticGraphGenerator for BarabasiAlbertGenerator {
    fn generate(&mut self) -> crate::Graph {
        BarabasiAlbertGeneratorGenerate(self.inner.pin_mut()).into()
    }
}