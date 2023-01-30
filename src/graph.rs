use cxx::UniquePtr;
use miette::{bail, IntoDiagnostic, Result};

use crate::bridge;

pub struct Graph {
    inner: UniquePtr<bridge::Graph>,
}

impl Graph {
    pub fn new(n: u64,
               weighted: bool,
               directed: bool,
               edges_indexed: bool) -> Self {
        Self {
            inner: bridge::NewGraph(n, weighted, directed, edges_indexed)
        }
    }

    pub fn add_edge(
        &mut self,
        u: u64,
        v: u64,
        ew: f64,
        check_multi_edge: bool,
    ) -> bool {
        self.inner.pin_mut().addEdge(u, v, ew, check_multi_edge)
    }

    pub fn add_node(&mut self) -> u64 {
        self.inner.pin_mut().addNode()
    }

    pub fn add_nodes(&mut self, number_of_new_nodes: u64) -> u64 {
        self.inner.pin_mut().addNodes(number_of_new_nodes)
    }

    pub fn check_consistency(&self) -> bool {
        self.inner.checkConsistency()
    }

    pub fn compact_edges(&mut self) {
        self.inner.pin_mut().compactEdges()
    }

    pub unsafe fn degree_unchecked(&self, v: u64) -> u64 {
        self.inner.degree(v)
    }

    pub fn degree(&self, v: u64) -> Result<u64> {
        if self.inner.hasNode(v) {
            Ok(unsafe {self.inner.degree(v)})
        } else {
            bail!("Node {} doesn't exist", v)
        }
    }

    pub unsafe fn degree_in_unchecked(&self, v: u64) -> u64 {
        self.inner.degreeIn(v)
    }

    pub fn degree_in(&self, v: u64) -> Result<u64> {
        if self.inner.hasNode(v) {
            Ok(unsafe {self.inner.degreeIn(v)})
        } else {
            bail!("Node {} doesn't exist", v)
        }
    }

    pub unsafe fn degree_out_unchecked(&self, v: u64) -> u64 {
        self.inner.degreeOut(v)
    }

    pub fn degree_out(&self, v: u64) -> Result<u64> {
        if self.inner.hasNode(v) {
            Ok(unsafe {self.inner.degreeOut(v)})
        } else {
            bail!("Node {} doesn't exist", v)
        }
    }

    pub fn edge_id(&self, u: u64, v: u64) -> Result<u64> {
        unsafe {self.inner.edgeId(u, v).into_diagnostic()}
    }

    pub fn has_edge(&self, u: u64, v: u64) -> bool {
        self.inner.hasEdge(u, v)
    }

    pub fn has_edge_ids(&self) -> bool {
        self.inner.hasEdgeIds()
    }

    pub fn has_node(&self, v: u64) -> bool {
        self.inner.hasNode(v)
    }

    pub unsafe fn increase_weight_unchecked(&mut self, u: u64, v: u64, ew: f64) -> Result<()> {
        self.inner.pin_mut().increaseWeight(u, v, ew).into_diagnostic()
    }

    pub fn increase_weight(&mut self, u: u64, v: u64, ew: f64) -> Result<()> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        if !self.inner.hasNode(v) {
            bail!("Node {} doesn't exist", v)
        }
        unsafe {self.inner.pin_mut().increaseWeight(u, v, ew).into_diagnostic()}
    }

    pub fn index_edges(&mut self, force: bool) {
        self.inner.pin_mut().indexEdges(force)
    }

    pub fn is_directed(&self) -> bool {
        self.inner.isDirected()
    }

    pub unsafe fn is_isolated_unchecked(&self, u: u64) -> Result<bool> {
        self.inner.isIsolated(u).into_diagnostic()
    }

    pub fn is_isolated(&self, u: u64) -> Result<bool> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        self.inner.isIsolated(u).into_diagnostic()
    }

    pub fn is_weighted(&self) -> bool {
        self.inner.isWeighted()
    }
}