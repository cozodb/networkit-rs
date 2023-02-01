use cxx::{CxxVector, UniquePtr};

use crate::bridge;
use crate::graph::Graph;
use miette::{bail, IntoDiagnostic, Result};

pub fn append(g: &mut Graph, g1: &Graph) {
    bridge::append(g.inner.pin_mut(), &g1.inner);
}

pub fn augment_graph(g: &mut Graph) -> u64 {
    bridge::augmentGraph(g.inner.pin_mut())
}

pub fn copy_nodes(g: &Graph) -> Graph {
    Graph {
        inner: bridge::GTCopyNodes(&g.inner),
    }
}

pub fn create_augmented_graph(g: &Graph) -> (Graph, u64) {
    let mut root = 0;
    let inner = bridge::GTCreateAugmentedGraph(g, &mut root);
    (Graph { inner }, root)
}

pub fn density(g: &Graph) -> f64 {
    bridge::density(&g.inner)
}

pub fn get_compacted_graph(g: &Graph, random: bool) -> Graph {
    let inner = bridge::GTGetCompactedGraph(g, random);
    Graph { inner }
}

pub fn volume(g: &Graph, nodes: &[u64]) -> f64 {
    bridge::GTVolume(g, nodes)
}

pub fn in_volume(g: &Graph, nodes: &[u64]) -> f64 {
    bridge::GTInVolume(g, nodes)
}

pub fn max_degree(g: &Graph) -> u64 {
    bridge::maxDegree(g)
}

pub fn max_in_degree(g: &Graph) -> u64 {
    bridge::maxInDegree(g)
}

pub fn max_weighted_degree(g: &Graph) -> f64 {
    bridge::maxWeightedDegree(g)
}

pub fn max_weighted_in_degree(g: &Graph) -> f64 {
    bridge::maxWeightedInDegree(g)
}

pub fn merge(g: &mut Graph, g1: &Graph) {
    bridge::merge(g.inner.pin_mut(), g1)
}

pub fn random_edge(g: &Graph, uniform: bool) -> Option<(u64, u64)> {
    let mut src = u64::MAX;
    let mut dst = u64::MAX;
    bridge::GTRandomEdge(g, uniform, &mut src, &mut dst);
    if src == u64::MAX || dst == u64::MAX {
        None
    } else {
        Some((src, dst))
    }
}

pub fn random_edges(g: &Graph, n: u64) -> impl Iterator<Item = (u64, u64)> {
    let mut src = vec![];
    let mut dst = vec![];
    bridge::GTRandomEdges(g, n, &mut src, &mut dst);
    src.into_iter().zip(dst.into_iter())
}
pub unsafe fn random_neighbour_unchecked(g: &Graph, u: u64) -> Option<u64> {
    match bridge::randomNeighbor(g, u) {
        u64::MAX => None,
        n => Some(n),
    }
}

pub fn random_neighbour(g: &Graph, u: u64) -> Result<Option<u64>> {
    if u >= g.upperEdgeIdBound() {
        bail!("Note {} out of bound", u)
    }
    Ok(unsafe { random_neighbour_unchecked(g, u) })
}

pub fn random_node(g: &Graph) -> Option<u64> {
    match bridge::randomNode(g) {
        u64::MAX => None,
        n => Some(n),
    }
}

pub fn random_nodes(g: &Graph, n: u64) -> impl Iterator<Item = u64> {
    let nodes = bridge::GTRandomNodes(g, n);
    NodeIter { nodes, at: 0 }
}

pub fn remove_edges_from_isolated_set(g: &mut Graph, nodes: &[u64]) {
    bridge::GTRemoveEdgesFromIsolatedSet(g.inner.pin_mut(), nodes)
}

pub fn size(g: &Graph) -> (u64, u64) {
    let mut n_size = 0;
    let mut e_size = 0;
    bridge::GTSize(g, &mut n_size, &mut e_size);
    (n_size, e_size)
}

pub fn sort_edges_by_weight(g: &mut Graph, descending: bool) {
    bridge::sortEdgesByWeight(g.inner.pin_mut(), descending)
}

pub fn subgraph_and_neighbours_from_nodes(
    g: &Graph,
    nodes: &[u64],
    include_out_neighbours: bool,
    include_in_neighbours: bool,
) -> Graph {
    Graph {
        inner: bridge::GTSubgraphAndNeighborsFromNodes(
            g,
            nodes,
            include_out_neighbours,
            include_in_neighbours,
        ),
    }
}

pub fn subgraph_from_nodes(g: &Graph, nodes: &[u64]) -> Graph {
    Graph {
        inner: bridge::GTSubgraphFromNodes(g, nodes),
    }
}

pub fn to_undirected(g: &Graph) -> Graph {
    Graph {
        inner: bridge::GTToUndirected(g),
    }
}

pub fn to_unweighted(g: &Graph) -> Graph {
    Graph {
        inner: bridge::GTToUnweighted(g),
    }
}

pub fn to_weighted(g: &Graph) -> Graph {
    Graph {
        inner: bridge::GTToWeighted(g),
    }
}

pub fn topological_sort(g: &Graph) -> impl Iterator<Item = u64> {
    NodeIter {
        nodes: bridge::GTTopologicalSort(g),
        at: 0,
    }
}

pub fn transpose(g: &Graph) -> Result<Graph> {
    let inner = bridge::GTTranspose(g).into_diagnostic()?;
    Ok(Graph { inner })
}

struct NodeIter {
    nodes: UniquePtr<CxxVector<u64>>,
    at: usize,
}

impl Iterator for NodeIter {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        match self.nodes.get(self.at) {
            None => None,
            Some(n) => {
                self.at += 1;
                Some(*n)
            }
        }
    }
}
