use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
    community::ValueIter,
    tools::NodeIter,
};

pub trait Centrality: Algorithm {
    fn centralization(&mut self) -> f64;
    fn maximum(&mut self) -> f64;
    fn ranking(&mut self) -> RankIter;
    fn score(&mut self, v: u64) -> f64;
    fn scores(&mut self) -> ValueIter;
}

pub struct RankIter {
    pub(crate) ks: Vec<u64>,
    pub(crate) vs: Vec<f64>,
    pub(crate) at: usize,
}

impl Iterator for RankIter {
    type Item = (u64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if self.at >= self.ks.len() {
            None
        } else {
            let k = self.ks[self.at];
            let v = self.vs[self.at];
            self.at += 1;
            Some((k, v))
        }
    }
}

pub struct ApproxBetweenness {
    inner: UniquePtr<bridge::ApproxBetweenness>,
}

impl ApproxBetweenness {
    pub fn new(
        g: &Graph,
        epsilon: Option<f64>,
        delta: Option<f64>,
        universal_constant: Option<f64>,
    ) -> Self {
        Self {
            inner: NewApproxBetweenness(
                g,
                epsilon.unwrap_or(0.01),
                delta.unwrap_or(0.1),
                universal_constant.unwrap_or(1.),
            ),
        }
    }
}

impl Centrality for ApproxBetweenness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproxBetweennessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxBetweennessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ApproxBetweenness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ApproxCloseness {
    inner: UniquePtr<bridge::ApproxCloseness>,
}

#[derive(Default)]
#[repr(u8)]
pub enum ApproxClosenessType {
    #[default]
    OutBound = 0,
    InBound = 1,
    Sum = 2,
}

impl ApproxCloseness {
    pub fn new(
        g: &Graph,
        n_samples: u64,
        epsilon: Option<f64>,
        normalized: bool,
        t: ApproxClosenessType,
    ) -> Self {
        Self {
            inner: NewApproxCloseness(g, n_samples, epsilon.unwrap_or(0.1), normalized, t as u8),
        }
    }

    pub fn get_square_error_estimates(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxClosenessGetSquareErrorEstimates(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Centrality for ApproxCloseness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproxClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxClosenessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ApproxCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct ApproxElectricalCloseness {
    inner: UniquePtr<bridge::ApproxElectricalCloseness>,
}

impl ApproxElectricalCloseness {
    pub fn new(g: &Graph, epsilon: Option<f64>, kappa: Option<f64>) -> Self {
        Self {
            inner: NewApproxElectricalCloseness(g, epsilon.unwrap_or(0.1), kappa.unwrap_or(0.3)),
        }
    }

    pub fn compute_exact_diagonal(&self, tol: Option<f64>) -> ValueIter {
        ValueIter {
            inner: ApproxElectricalClosenessComputeExactDiagonal(&self.inner, tol.unwrap_or(1e-9)),
            at: 0,
        }
    }
    pub fn get_diagonal(&self) -> ValueIter {
        ValueIter {
            inner: ApproxElectricalClosenessGetDiagonal(&self.inner),
            at: 0,
        }
    }
}

impl Centrality for ApproxElectricalCloseness {
    fn centralization(&mut self) -> f64 {
        self.inner.pin_mut().centralization()
    }

    fn maximum(&mut self) -> f64 {
        self.inner.pin_mut().maximum()
    }

    fn ranking(&mut self) -> RankIter {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproxElectricalClosenessRanking(self.inner.pin_mut(), &mut ks, &mut vs);
        RankIter { ks, vs, at: 0 }
    }

    fn score(&mut self, node: u64) -> f64 {
        self.inner.pin_mut().score(node)
    }

    fn scores(&mut self) -> ValueIter {
        ValueIter {
            inner: ApproxElectricalClosenessScores(self.inner.pin_mut()),
            at: 0,
        }
    }
}

impl Algorithm for ApproxElectricalCloseness {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}
