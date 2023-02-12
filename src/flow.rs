use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
    tools::NodeIter,
};

pub struct EdmondsKarp {
    inner: UniquePtr<bridge::EdmondsKarp>,
}

impl EdmondsKarp {
    pub fn new(g: &crate::Graph, source: u64, sink: u64) -> Self {
        Self {
            inner: NewEdmondsKarp(g, source, sink),
        }
    }
    pub fn get_max_flow(&self) -> f64 {
        self.inner.getMaxFlow()
    }

    pub fn get_source_set(&self) -> NodeIter {
        NodeIter {
            at: 0,
            nodes: EdmondsKarpGetSourceSet(&self.inner),
        }
    }
    pub fn get_flow(&self, u: u64, v: u64) -> f64 {
        self.inner.getFlow(u, v)
    }
    pub fn get_edge_flow(&self, e: u64) -> f64 {
        self.inner.getEdgeFlow(e)
    }
    pub fn get_flow_vector(&self) -> Vec<f64> {
        EdmondsKarpGetFlowVector(&self.inner)
            .iter()
            .cloned()
            .collect()
    }
}

impl Algorithm for EdmondsKarp {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}
