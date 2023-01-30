use autocxx::prelude::*;

// use all the main autocxx functions

include_cpp! {
    #include "networkit/centrality/Centrality.hpp"
    #include "networkit/centrality/PageRank.hpp"
    #include "extra.h"

    safety!(unsafe)
    generate!("NetworKit::Graph")
    generate!("NetworKit::Algorithm")
    generate!("NetworKit::Centrality")
    generate!("NetworKit::PageRank")
    generate!("NetworKit::PageRank_Norm")
    generate!("NetworKit::PageRank_SinkHandling")
    generate!("NetworKit::NewPageRank")
}

fn main() {
    use crate::ffi::NetworKit::*;

    println!("hello, world!!!");
    let g = Graph::new(10, false, false, false).within_unique_ptr();
    println!("{}", g.numberOfNodes());
}
