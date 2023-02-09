use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
};



pub struct Assortativity {
    inner: UniquePtr<bridge::Assortativity>,
}

impl Assortativity {
    pub fn new(g: &crate::Graph, attributes: &[f64]) -> Self {
        Self {
            inner: NewAssortativity(g, attributes),
        }
    }
    pub fn get_coefficient(&self) -> f64 {
        self.inner.getCoefficient()
    }
}

impl Algorithm for Assortativity {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

