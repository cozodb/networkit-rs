use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
};

pub struct MaximalCliques {
    inner: UniquePtr<bridge::MaximalCliques>,
}

impl MaximalCliques {
    pub fn new(g: &crate::Graph, maximum_only: bool) -> Self {
        Self {
            inner: NewMaximalCliques(g, maximum_only),
        }
    }
    pub fn get_cliques(&mut self) -> impl Iterator<Item = (u64, u64)> {
        let mut cliques = vec![];
        let mut nodes = vec![];
        MaximalCliquesGetCliques(self.inner.pin_mut(), &mut cliques, &mut nodes);
        cliques.into_iter().zip(nodes.into_iter())
    }
}

impl Algorithm for MaximalCliques {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}
