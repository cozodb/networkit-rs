use cxx::UniquePtr;
use miette::{bail, IntoDiagnostic, Result};

use crate::bridge::{self, *};

pub struct Graph {
    inner: UniquePtr<bridge::Graph>,
}

impl Graph {
    pub fn new(n: u64, weighted: bool, directed: bool, edges_indexed: bool) -> Self {
        Self {
            inner: bridge::NewGraph(n, weighted, directed, edges_indexed),
        }
    }

    pub fn add_edge(&mut self, u: u64, v: u64, ew: f64, check_multi_edge: bool) -> bool {
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
            Ok(unsafe { self.inner.degree(v) })
        } else {
            bail!("Node {} doesn't exist", v)
        }
    }

    pub unsafe fn degree_in_unchecked(&self, v: u64) -> u64 {
        self.inner.degreeIn(v)
    }

    pub fn degree_in(&self, v: u64) -> Result<u64> {
        if self.inner.hasNode(v) {
            Ok(unsafe { self.inner.degreeIn(v) })
        } else {
            bail!("Node {} doesn't exist", v)
        }
    }

    pub unsafe fn degree_out_unchecked(&self, v: u64) -> u64 {
        self.inner.degreeOut(v)
    }

    pub fn degree_out(&self, v: u64) -> Result<u64> {
        if self.inner.hasNode(v) {
            Ok(unsafe { self.inner.degreeOut(v) })
        } else {
            bail!("Node {} doesn't exist", v)
        }
    }

    pub fn edge_id(&self, u: u64, v: u64) -> Result<u64> {
        unsafe { self.inner.edgeId(u, v).into_diagnostic() }
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
        self.inner
            .pin_mut()
            .increaseWeight(u, v, ew)
            .into_diagnostic()
    }

    pub fn increase_weight(&mut self, u: u64, v: u64, ew: f64) -> Result<()> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        if !self.inner.hasNode(v) {
            bail!("Node {} doesn't exist", v)
        }
        unsafe {
            self.inner
                .pin_mut()
                .increaseWeight(u, v, ew)
                .into_diagnostic()
        }
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

    pub fn number_of_edges(&self) -> u64 {
        self.inner.numberOfEdges()
    }

    pub fn number_of_nodes(&self) -> u64 {
        self.inner.numberOfNodes()
    }

    pub fn number_of_self_loops(&self) -> u64 {
        self.inner.numberOfSelfLoops()
    }

    pub fn remove_all_edges(&mut self) {
        self.inner.pin_mut().removeAllEdges()
    }

    pub fn remove_edge(&mut self, u: u64, v: u64) -> Result<()> {
        self.inner.pin_mut().removeEdge(u, v).into_diagnostic()
    }

    pub fn remove_multi_edges(&mut self) {
        self.inner.pin_mut().removeMultiEdges()
    }

    pub unsafe fn remove_node_unchecked(&mut self, u: u64) {
        self.inner.pin_mut().removeNode(u)
    }

    pub fn remove_node(&mut self, u: u64) -> Result<()> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        unsafe { self.inner.pin_mut().removeNode(u) };
        Ok(())
    }

    pub fn remove_self_loops(&mut self) {
        self.inner.pin_mut().removeSelfLoops()
    }

    pub unsafe fn restore_node_unchecked(&mut self, u: u64) {
        self.inner.pin_mut().restoreNode(u)
    }

    pub fn restore_node(&mut self, u: u64) -> Result<()> {
        if u >= self.inner.upperNodeIdBound() {
            bail!("Node {} out of bound", u)
        }
        unsafe { self.inner.pin_mut().restoreNode(u) }
        Ok(())
    }
    pub unsafe fn set_weight_unchecked(&mut self, u: u64, v: u64, ew: f64) -> Result<()> {
        self.inner.pin_mut().setWeight(u, v, ew).into_diagnostic()
    }

    pub fn set_weight(&mut self, u: u64, v: u64, ew: f64) -> Result<()> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        if !self.inner.hasNode(v) {
            bail!("Node {} doesn't exist", v)
        }
        unsafe { self.inner.pin_mut().setWeight(u, v, ew).into_diagnostic() }
    }

    pub fn sort_edges(&mut self) {
        self.inner.pin_mut().sortEdges()
    }

    pub unsafe fn swap_edge_unchecked(&mut self, s1: u64, t1: u64, s2: u64, t2: u64) {
        self.inner.pin_mut().swapEdge(s1, t1, s2, t2)
    }

    pub fn swap_edge(&mut self, s1: u64, t1: u64, s2: u64, t2: u64) -> Result<()> {
        for u in [s1, t1, s2, t2] {
            if !self.inner.hasNode(u) {
                bail!("Node {} doesn't exist", u)
            }
        }
        unsafe { self.inner.pin_mut().swapEdge(s1, t1, s2, t2) }
        Ok(())
    }
    pub fn total_edge_weight(&self) -> f64 {
        self.inner.totalEdgeWeight()
    }
    pub fn upper_edge_id_bound(&self) -> u64 {
        self.inner.upperEdgeIdBound()
    }
    pub fn upper_node_id_bound(&self) -> u64 {
        self.inner.upperNodeIdBound()
    }

    pub unsafe fn weight_unchecked(&self, u: u64, v: u64) -> f64 {
        self.inner.weight(u, v)
    }

    pub fn weight(&self, u: u64, v: u64) -> Result<f64> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        if !self.inner.hasNode(v) {
            bail!("Node {} doesn't exist", v)
        }
        Ok(unsafe { self.inner.weight(u, v) })
    }

    pub unsafe fn weighted_degree_unchecked(
        self: &Graph,
        u: u64,
        count_self_loops_twice: bool,
    ) -> f64 {
        self.inner.weightedDegree(u, count_self_loops_twice)
    }

    pub fn weighted_degree(&self, u: u64, count_self_loops_twice: bool) -> Result<f64> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        Ok(unsafe { self.inner.weightedDegree(u, count_self_loops_twice) })
    }
    pub unsafe fn weighted_degree_in_unchecked(
        self: &Graph,
        u: u64,
        count_self_loops_twice: bool,
    ) -> f64 {
        self.inner.weightedDegreeIn(u, count_self_loops_twice)
    }

    pub fn weighted_degree_in(&self, u: u64, count_self_loops_twice: bool) -> Result<f64> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }
        Ok(unsafe { self.inner.weightedDegreeIn(u, count_self_loops_twice) })
    }

    pub fn iter_nodes<'a>(&'a self) -> impl Iterator<Item = u64> + 'a {
        struct It(UniquePtr<GraphNodeIter>);

        impl Iterator for It {
            type Item = u64;
            fn next(&mut self) -> Option<Self::Item> {
                let mut cur = 0;
                if self.0.pin_mut().advance(&mut cur) {
                    Some(cur)
                } else {
                    None
                }
            }
        }

        It(NewGraphNodeIter(&self.inner))
    }
    pub fn iter_edges<'a>(&'a self) -> impl Iterator<Item = (u64, u64)> + 'a {
        struct It(UniquePtr<GraphEdgeIter>);

        impl Iterator for It {
            type Item = (u64, u64);
            fn next(&mut self) -> Option<Self::Item> {
                let mut src = 0;
                let mut dst = 0;
                if self.0.pin_mut().advance(&mut src, &mut dst) {
                    Some((src, dst))
                } else {
                    None
                }
            }
        }

        It(NewGraphEdgeIter(&self.inner))
    }

    pub fn iter_edges_weight<'a>(&'a self) -> impl Iterator<Item = (u64, u64, f64)> + 'a {
        struct It(UniquePtr<GraphEdgeWeightIter>);

        impl Iterator for It {
            type Item = (u64, u64, f64);
            fn next(&mut self) -> Option<Self::Item> {
                let mut src = 0;
                let mut dst = 0;
                let mut wt = 0.0;
                if self.0.pin_mut().advance(&mut src, &mut dst, &mut wt) {
                    Some((src, dst, wt))
                } else {
                    None
                }
            }
        }

        It(NewGraphEdgeWeightIter(&self.inner))
    }
    pub fn iter_neighbours<'a>(&'a self, u: u64) -> Result<impl Iterator<Item = u64> + 'a> {
        self.iter_neighbours_impl(u, false)
    }
    pub fn iter_in_neighbours<'a>(&'a self, u: u64) -> Result<impl Iterator<Item = u64> + 'a> {
        self.iter_neighbours_impl(u, true)
    }
    fn iter_neighbours_impl<'a>(
        &'a self,
        u: u64,
        in_neighbours: bool,
    ) -> Result<impl Iterator<Item = u64> + 'a> {
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }

        struct It(UniquePtr<GraphNeighbourIter>);
        impl Iterator for It {
            type Item = u64;
            fn next(&mut self) -> Option<Self::Item> {
                let mut cur = 0;
                if self.0.pin_mut().advance(&mut cur) {
                    Some(cur)
                } else {
                    None
                }
            }
        }
        Ok(It(unsafe {
            NewGraphNeighbourIter(&self.inner, u, in_neighbours)
        }))
    }
    pub fn iter_neighbours_weight<'a>(
        &'a self,
        u: u64,
    ) -> Result<impl Iterator<Item = (u64, f64)> + 'a> {
        self.iter_neighbours_weight_impl(u, false)
    }
    pub fn iter_in_neighbours_weight<'a>(
        &'a self,
        u: u64,
    ) -> Result<impl Iterator<Item = (u64, f64)> + 'a> {
        self.iter_neighbours_weight_impl(u, true)
    }
    fn iter_neighbours_weight_impl<'a>(
        &'a self,
        u: u64,
        in_neighbours: bool,
    ) -> Result<impl Iterator<Item = (u64, f64)> + 'a> {
        if !self.inner.isWeighted() {
            bail!("Graph is unweighted")
        }
        if !self.inner.hasNode(u) {
            bail!("Node {} doesn't exist", u)
        }

        struct It(UniquePtr<GraphNeighbourWeightIter>);
        impl Iterator for It {
            type Item = (u64, f64);
            fn next(&mut self) -> Option<Self::Item> {
                let mut cur = 0;
                let mut wt = 0.0;
                if self.0.pin_mut().advance(&mut cur, &mut wt) {
                    Some((cur, wt))
                } else {
                    None
                }
            }
        }
        Ok(It(unsafe {
            NewGraphNeighbourWeightIter(&self.inner, u, in_neighbours).into_diagnostic()?
        }))
    }
}
