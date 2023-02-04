use networkit_rs::base::Algorithm;
use networkit_rs::community::{community_graph, CommunityDetector, PLM};
use networkit_rs::graph::Graph;
use networkit_rs::tools::to_undirected;

#[test]
fn creation() {
    let mut g = Graph::new(10, false, true, false);
    g.add_edge(0, 2, None, false);
    g.add_edge(1, 3, None, false);
    g.add_edge(9, 0, None, false);
    println!("Nodes");
    for u in g.iter_nodes() {
        println!("{u}");
    }
    println!("Edges");
    for (u, v, wt) in g.iter_edges_weight() {
        println!("{u} -> {v}: {wt}")
    }
    println!("Neigbours");
    for u in g.iter_neighbours(0).unwrap() {
        println!("{u}")
    }
    println!("InNeigbours");
    for u in g.iter_in_neighbours(0).unwrap() {
        println!("{u}")
    }

    let g = to_undirected(&g);

    let mut plm = PLM::new(&g, None, None, None, None, None, None).unwrap();
    plm.run().unwrap();
    let p = plm.get_partition();
    let ng = community_graph(&g, &p).unwrap();

    println!("Community edges");
    for (u, v, wt) in ng.iter_edges_weight() {
        println!("{u} -> {v}: {wt}")
    }
}
