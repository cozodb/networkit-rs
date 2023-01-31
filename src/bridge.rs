extern crate openmp_sys;

pub(crate) use ffi::*;

//
// // use all the main autocxx functions
//
// include_cpp! {
//     #include "networkit/centrality/Centrality.hpp"
//     #include "networkit/centrality/PageRank.hpp"
//     #include "extra.h"
//
//     safety!(unsafe)
//     generate!("NetworKit::Graph")
//     generate!("NetworKit::Algorithm")
//     generate!("NetworKit::Centrality")
//     generate!("NetworKit::PageRank")
//     generate!("NetworKit::PageRank_Norm")
//     generate!("NetworKit::PageRank_SinkHandling")
//     generate!("NetworKit::NewPageRank")
// }

#[cxx::bridge(namespace = "NetworKit")]
mod ffi {
    unsafe extern "C++" {
        include!("bridge.h");

        // GRAPH

        pub type Graph;

        pub fn NewGraph(
            n: u64,
            weighted: bool,
            directed: bool,
            edges_indexed: bool,
        ) -> UniquePtr<Graph>;
        pub fn addEdge(
            self: Pin<&mut Graph>,
            u: u64,
            v: u64,
            ew: f64,
            check_multi_edge: bool,
        ) -> bool;
        fn addNode(self: Pin<&mut Graph>) -> u64;
        fn addNodes(self: Pin<&mut Graph>, number_of_new_nodes: u64) -> u64;
        fn checkConsistency(self: &Graph) -> bool;
        fn compactEdges(self: Pin<&mut Graph>);
        unsafe fn degree(self: &Graph, v: u64) -> u64;
        unsafe fn degreeIn(self: &Graph, v: u64) -> u64;
        unsafe fn degreeOut(self: &Graph, v: u64) -> u64;
        unsafe fn edgeId(self: &Graph, u: u64, v: u64) -> Result<u64>;
        fn hasEdge(self: &Graph, u: u64, v: u64) -> bool;
        fn hasEdgeIds(self: &Graph) -> bool;
        fn hasNode(self: &Graph, v: u64) -> bool;
        unsafe fn increaseWeight(self: Pin<&mut Graph>, u: u64, v: u64, ew: f64) -> Result<()>;
        fn indexEdges(self: Pin<&mut Graph>, force: bool);
        fn isDirected(self: &Graph) -> bool;
        fn isIsolated(self: &Graph, u: u64) -> Result<bool>;
        fn isWeighted(self: &Graph) -> bool;

        fn numberOfEdges(self: &Graph) -> u64;
        fn numberOfNodes(self: &Graph) -> u64;
        fn numberOfSelfLoops(self: &Graph) -> u64;

        fn removeAllEdges(self: Pin<&mut Graph>);
        fn removeEdge(self: Pin<&mut Graph>, u: u64, v: u64) -> Result<()>;
        fn removeMultiEdges(self: Pin<&mut Graph>);
        unsafe fn removeNode(self: Pin<&mut Graph>, u: u64);
        fn removeSelfLoops(self: Pin<&mut Graph>);
        unsafe fn restoreNode(self: Pin<&mut Graph>, u: u64);
        unsafe fn setWeight(self: Pin<&mut Graph>, u: u64, v: u64, ew: f64) -> Result<()>;

        fn sortEdges(self: Pin<&mut Graph>);
        unsafe fn swapEdge(self: Pin<&mut Graph>, s1: u64, t1: u64, s2: u64, t2: u64);

        fn totalEdgeWeight(self: &Graph) -> f64;

        fn upperEdgeIdBound(self: &Graph) -> u64;
        fn upperNodeIdBound(self: &Graph) -> u64;

        unsafe fn weight(self: &Graph, u: u64, v: u64) -> f64;

        unsafe fn weightedDegree(self: &Graph, u: u64, count_self_loops_twice: bool) -> f64;
        unsafe fn weightedDegreeIn(self: &Graph, u: u64, count_self_loops_twice: bool) -> f64;

        pub type GraphNodeIter;

        fn NewGraphNodeIter(g: &Graph) -> UniquePtr<GraphNodeIter>;
        fn advance(self: Pin<&mut GraphNodeIter>, u: &mut u64) -> bool;

        pub type GraphEdgeIter;
        fn NewGraphEdgeIter(g: &Graph) -> UniquePtr<GraphEdgeIter>;
        fn advance(self: Pin<&mut GraphEdgeIter>, u: &mut u64, v: &mut u64) -> bool;

        pub type GraphEdgeWeightIter;
        fn NewGraphEdgeWeightIter(g: &Graph) -> UniquePtr<GraphEdgeWeightIter>;
        fn advance(
            self: Pin<&mut GraphEdgeWeightIter>,
            u: &mut u64,
            v: &mut u64,
            wt: &mut f64,
        ) -> bool;

        pub type GraphNeighbourIter;
        unsafe fn NewGraphNeighbourIter(
            g: &Graph,
            u: u64,
            in_neighbours: bool,
        ) -> UniquePtr<GraphNeighbourIter>;
        fn advance(self: Pin<&mut GraphNeighbourIter>, u: &mut u64) -> bool;

        pub type GraphNeighbourWeightIter;
        unsafe fn NewGraphNeighbourWeightIter(
            g: &Graph,
            u: u64,
            in_neighbours: bool,
        ) -> Result<UniquePtr<GraphNeighbourWeightIter>>;
        fn advance(self: Pin<&mut GraphNeighbourWeightIter>, u: &mut u64, wt: &mut f64) -> bool;
    }
}
