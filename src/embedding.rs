use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
};

pub struct Node2Vec {
    inner: UniquePtr<bridge::Node2Vec>,
    n_nodes: usize,
    dim: usize,
}

#[derive(Clone)]
pub struct NodeFeatures {
    data: Vec<f64>,
    n_nodes: usize,
    dim: usize,
}

impl NodeFeatures {
    pub fn iter<'a>(&'a self) -> impl Iterator<Item = &'a [f64]> {
        (0..self.n_nodes).map(|idx| self.get(idx))
    }
    pub fn get(&self, idx: usize) -> &[f64] {
        &self.data[idx * self.dim..idx * self.dim + self.dim]
    }
}

impl Node2Vec {
    pub fn new(
        g: &crate::Graph,
        p: Option<f64>,
        q: Option<f64>,
        l: Option<u64>,
        n: Option<u64>,
        d: Option<u64>,
    ) -> Self {
        let d = d.unwrap_or(128);
        let n_nodes = g.number_of_nodes() as usize;
        Self {
            inner: NewNode2Vec(
                g,
                p.unwrap_or(1.),
                q.unwrap_or(1.),
                l.unwrap_or(80),
                n.unwrap_or(10),
                d,
            ),
            n_nodes,
            dim: d as usize,
        }
    }
    pub fn get_features(&self) -> NodeFeatures {
        let mut data = Vec::with_capacity(self.n_nodes * self.dim);
        Node2VecGetFeatures(&self.inner, &mut data);
        NodeFeatures {
            data,
            n_nodes: self.n_nodes,
            dim: self.dim,
        }
    }
}

impl Algorithm for Node2Vec {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}
