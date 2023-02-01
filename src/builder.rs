use cxx::UniquePtr;

use crate::bridge::{self, *};
use crate::graph;
use miette::{bail, Result};

pub struct GraphBuilder {
    inner: UniquePtr<bridge::GraphBuilder>,
}

unsafe impl Send for GraphBuilder {}

impl GraphBuilder {
    pub fn new(n: u64, weighted: bool, directed: bool) -> Self {
        Self {
            inner: NewGraphBuilder(n, weighted, directed),
        }
    }
    pub fn complete_graph(&mut self, parallel: bool) -> graph::Graph {
        let g = GraphBuilderCompleteGraph(self.inner.pin_mut(), parallel);
        graph::Graph { inner: g }
    }
    pub fn reset(&mut self, n: u64) {
        self.inner.pin_mut().reset(n);
    }
    pub fn is_weighted(&self) -> bool {
        self.inner.isWeighted()
    }
    pub fn is_directed(&self) -> bool {
        self.inner.isDirected()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.isEmpty()
    }
    pub fn number_of_nodes(&self) -> u64 {
        self.inner.numberOfNodes()
    }
    pub fn upper_node_id_bound(&self) -> u64 {
        self.inner.upperNodeIdBound()
    }
    pub fn add_node(&mut self) -> u64 {
        self.inner.pin_mut().addNode()
    }
    pub unsafe fn add_half_edge_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().addHalfEdge(u, v, ew.unwrap_or(1.))
    }
    pub fn add_half_edge(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.add_half_edge_unchecked(u, v, ew) })
    }
    pub unsafe fn add_half_out_edge_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().addHalfOutEdge(u, v, ew.unwrap_or(1.))
    }
    pub fn add_half_out_edge(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.add_half_out_edge_unchecked(u, v, ew) })
    }
    pub unsafe fn add_half_in_edge_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().addHalfInEdge(u, v, ew.unwrap_or(1.))
    }
    pub fn add_half_in_edge(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.add_half_in_edge_unchecked(u, v, ew) })
    }
    pub unsafe fn set_weight_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().setWeight(u, v, ew.unwrap_or(1.))
    }
    pub fn set_weight(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.set_weight_unchecked(u, v, ew) })
    }
    pub unsafe fn set_out_weight_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().setOutWeight(u, v, ew.unwrap_or(1.))
    }
    pub fn set_out_weight(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.set_out_weight_unchecked(u, v, ew) })
    }
    pub unsafe fn set_in_weight_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().setInWeight(u, v, ew.unwrap_or(1.))
    }
    pub fn set_in_weight(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.set_in_weight_unchecked(u, v, ew) })
    }
    pub unsafe fn increase_weight_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner.pin_mut().increaseWeight(u, v, ew.unwrap_or(1.))
    }
    pub fn increase_weight(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.increase_weight_unchecked(u, v, ew) })
    }
    pub unsafe fn increase_out_weight_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner
            .pin_mut()
            .increaseOutWeight(u, v, ew.unwrap_or(1.))
    }
    pub fn increase_out_weight(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.increase_out_weight_unchecked(u, v, ew) })
    }
    pub unsafe fn increase_in_weight_unchecked(&mut self, u: u64, v: u64, ew: Option<f64>) {
        self.inner
            .pin_mut()
            .increaseInWeight(u, v, ew.unwrap_or(1.))
    }
    pub fn increase_in_weight(&mut self, u: u64, v: u64, ew: Option<f64>) -> Result<()> {
        for n in [u, v] {
            if n >= self.upper_node_id_bound() {
                bail!("Node out of bound: {}", n);
            }
        }
        Ok(unsafe { self.increase_in_weight_unchecked(u, v, ew) })
    }
}
