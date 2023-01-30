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

        pub type Graph;

        pub fn NewGraph(n: u64, weighted: bool, directed: bool, edges_indexed: bool) -> UniquePtr<Graph>;
        pub fn addEdge(
            self: Pin<&mut Graph>,
            u: u64,
            v: u64,
            ew: f64,
            checkMultiEdge: bool,
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

        // TODO iterators

        // pub type Graph_NodeRange;
        //
        // fn NewNodeRange(g: &Graph) -> UniquePtr<Graph_NodeRange>;
        //
        // pub type Graph_EdgeRange;
        //
        // fn NewEdgeRange(g: &Graph) -> UniquePtr<Graph_EdgeRange>;

        fn numberOfEdges(self: &Graph) -> u64;
        fn numberOfNodes(self: &Graph) -> u64;
        fn numberOfSelfLoops(self: &Graph) -> u64;

        fn removeAllEdges(self: Pin<&mut Graph>);
        fn removeEdge(self: Pin<&mut Graph>, u: u64, v: u64) -> Result<()>;
        fn removeMultiEdges(self: Pin<&mut Graph>);
        unsafe fn removeNode(self: Pin<&mut Graph>, u: u64);
        fn removeSelfLoops(self: Pin<&mut Graph>);

        // TODO

        fn upperNodeIdBound(self: &Graph) -> u64;
    }
}