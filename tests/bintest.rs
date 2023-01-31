use networkit_rs::graph::Graph;

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
}
