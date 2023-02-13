use cxx::{CxxVector, UniquePtr};
use miette::IntoDiagnostic;

use crate::{
    base::{Algorithm, GraphEvent, GraphEventType},
    bridge::{self, *},
    tools::NodeIter,
};

pub trait StaticGraphGenerator {
    fn generate(&mut self) -> crate::Graph;
}

pub trait DynamicGraphGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent>;
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

pub struct ChungLuGenerator {
    inner: UniquePtr<bridge::ChungLuGenerator>,
    _seq: UniquePtr<CxxVector<u64>>,
}

impl ChungLuGenerator {
    pub fn new(degree_sequence: &[u64]) -> Self {
        let seq = MakeCountVector(degree_sequence);
        let inner = NewChungLuGenerator(&seq);
        Self { inner, _seq: seq }
    }
}

impl StaticGraphGenerator for ChungLuGenerator {
    fn generate(&mut self) -> crate::Graph {
        ChungLuGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct ClusteredRandomGraphGenerator {
    inner: UniquePtr<bridge::ClusteredRandomGraphGenerator>,
}

impl ClusteredRandomGraphGenerator {
    pub fn new(n: u64, k: u64, p_intra: f64, p_inter: f64) -> Self {
        let inner = NewClusteredRandomGraphGenerator(n, k, p_intra, p_inter);
        Self { inner }
    }
    pub fn get_communities(&mut self) -> crate::Partition {
        ClusteredRandomGraphGeneratorGetCommunities(self.inner.pin_mut()).into()
    }
}

impl StaticGraphGenerator for ClusteredRandomGraphGenerator {
    fn generate(&mut self) -> crate::Graph {
        ClusteredRandomGraphGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct DorogovtsevMendesGenerator {
    inner: UniquePtr<bridge::DorogovtsevMendesGenerator>,
}

impl DorogovtsevMendesGenerator {
    pub fn new(n_nodes: u64) -> Self {
        let inner = NewDorogovtsevMendesGenerator(n_nodes);
        Self { inner }
    }
}

impl StaticGraphGenerator for DorogovtsevMendesGenerator {
    fn generate(&mut self) -> crate::Graph {
        DorogovtsevMendesGeneratorGenerate(self.inner.pin_mut()).into()
    }
}

pub struct DynamicDorogovtsevMendesGenerator {
    inner: UniquePtr<bridge::DynamicDorogovtsevMendesGenerator>,
}

impl DynamicDorogovtsevMendesGenerator {
    pub fn new() -> Self {
        let inner = NewDynamicDorogovtsevMendesGenerator();
        Self { inner }
    }
}

impl DynamicGraphGenerator for DynamicDorogovtsevMendesGenerator {
    fn generate(&mut self, n_steps: u64) -> Vec<crate::base::GraphEvent> {
        let mut typs = vec![];
        let mut us = vec![];
        let mut vs = vec![];
        let mut ws = vec![];
        DynamicDorogovtsevMendesGeneratorGenerate(
            self.inner.pin_mut(),
            n_steps,
            &mut typs,
            &mut us,
            &mut vs,
            &mut ws,
        );
        typs.into_iter()
            .zip(us)
            .zip(vs)
            .zip(ws)
            .map(|(((t, u), v), w)| GraphEvent {
                kind: GraphEventType::from(t),
                u,
                v,
                ew: w,
            })
            .collect()
    }
}
