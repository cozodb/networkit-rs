use networkit_rs::graph::Graph;

#[test]
fn creation() {
    let mut g = Graph::new(10, false, false, false);
    g.index_edges(false);
    let d = g.degree_in(100000000).unwrap();
    println!("{d}");
    let x = g.edge_id(1, 2).unwrap();
    println!("{x}");
}