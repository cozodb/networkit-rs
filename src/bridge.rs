extern crate openmp_sys;

pub(crate) use ffi::*;

#[cxx::bridge]
mod ffi {

    #[namespace = "NetworKit"]
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
        pub fn CopyGraph(g: &Graph) -> UniquePtr<Graph>;
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

        // GRAPH BUILDER

        pub type GraphBuilder;
        fn NewGraphBuilder(n: u64, weighted: bool, directed: bool) -> UniquePtr<GraphBuilder>;
        fn reset(self: Pin<&mut GraphBuilder>, n: u64);
        fn isWeighted(self: &GraphBuilder) -> bool;
        fn isDirected(self: &GraphBuilder) -> bool;
        fn isEmpty(self: &GraphBuilder) -> bool;
        fn numberOfNodes(self: &GraphBuilder) -> u64;
        fn upperNodeIdBound(self: &GraphBuilder) -> u64;
        fn addNode(self: Pin<&mut GraphBuilder>) -> u64;
        unsafe fn addHalfEdge(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn addHalfOutEdge(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn addHalfInEdge(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        // unsafe fn swapNeighborhood: not needed
        unsafe fn setWeight(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn setOutWeight(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn setInWeight(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn increaseWeight(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn increaseOutWeight(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        unsafe fn increaseInWeight(self: Pin<&mut GraphBuilder>, u: u64, v: u64, ew: f64);
        // completeGraph
        fn GraphBuilderCompleteGraph(
            builder: Pin<&mut GraphBuilder>,
            parallel: bool,
        ) -> UniquePtr<Graph>;
        // iterators for builders are omitted

        // PARTITION
        pub type Partition;
        pub fn NewPartition(z: u64) -> UniquePtr<Partition>;
        fn CopyPartition(p: &Partition) -> UniquePtr<Partition>;
        fn addToSubset(self: Pin<&mut Partition>, s: u64, e: u64);
        fn allToSingletons(self: Pin<&mut Partition>);
        fn compact(self: Pin<&mut Partition>, use_turbo: bool);
        fn contains(self: &Partition, e: u64) -> bool;
        fn extend(self: Pin<&mut Partition>) -> u64;
        fn PTGetMembers(p: &Partition, s: u64, rs: &mut Vec<u64>);
        fn PTGetName(p: &Partition) -> UniquePtr<CxxString>;
        fn PTGetSubsetIds(p: &Partition, rs: &mut Vec<u64>);
        fn getVector(self: &Partition) -> &CxxVector<u64>;
        fn inSameSubset(self: &Partition, e1: u64, e2: u64) -> bool;
        fn lowerBound(self: &Partition) -> u64;
        fn mergeSubsets(self: Pin<&mut Partition>, s: u64, t: u64) -> u64;
        fn moveToSubset(self: Pin<&mut Partition>, s: u64, e: u64);
        fn numberOfElements(self: &Partition) -> u64;
        fn numberOfSubsets(self: &Partition) -> u64;
        fn PTSetName(p: Pin<&mut Partition>, name: &str);
        fn setUpperBound(self: Pin<&mut Partition>, upper: u64);
        fn subsetOf(self: &Partition, e: u64) -> u64;
        fn PTSubsetSizeMap(p: &Partition, ks: &mut Vec<u64>, sz: &mut Vec<u64>);
        fn PTSubsetSizes(p: &Partition) -> UniquePtr<CxxVector<u64>>;
        fn toSingleton(self: Pin<&mut Partition>, e: u64);
        fn upperBound(self: &Partition) -> u64;

        // COVER
        pub type Cover;
        pub fn NewCover() -> UniquePtr<Cover>;
        pub fn NewCoverWithSize(z: u64) -> UniquePtr<Cover>;
        pub fn NewCoverFromPartition(p: &Partition) -> UniquePtr<Cover>;
        pub fn CopyCover(c: &Cover) -> UniquePtr<Cover>;
        fn addToSubset(self: Pin<&mut Cover>, s: u64, e: u64);
        fn allToSingletons(self: Pin<&mut Cover>);
        fn contains(self: &Cover, e: u64) -> bool;
        fn extend(self: Pin<&mut Cover>) -> u64;
        fn CVGetMembers(c: &Cover, s: u64, rs: &mut Vec<u64>);
        fn CVGetSubsetIds(c: &Cover, rs: &mut Vec<u64>);
        fn inSameSubset(self: &Cover, e1: u64, e2: u64) -> bool;
        fn lowerBound(self: &Cover) -> u64;
        fn mergeSubsets(self: Pin<&mut Cover>, s: u64, t: u64);
        fn moveToSubset(self: Pin<&mut Cover>, s: u64, e: u64);
        fn numberOfElements(self: &Cover) -> u64;
        fn numberOfSubsets(self: &Cover) -> u64;
        fn removeFromSubset(self: Pin<&mut Cover>, s: u64, e: u64);
        fn setUpperBound(self: Pin<&mut Cover>, upper: u64);
        fn CVSubsetSizeMap(c: &Cover, ks: &mut Vec<u64>, sz: &mut Vec<u64>);
        fn CVSubsetSizes(c: &Cover) -> UniquePtr<CxxVector<u64>>;
        fn CVSubsetsOf(c: &Cover, e: u64) -> UniquePtr<CxxVector<u64>>;
        fn toSingleton(self: Pin<&mut Cover>, e: u64) -> u64;
        fn upperBound(self: &Cover) -> u64;

        // COMMUNITY
        type AdjustedRandMeasure;
        pub fn NewAdjustedRandMeasure() -> UniquePtr<AdjustedRandMeasure>;
        pub fn getDissimilarity(
            self: Pin<&mut AdjustedRandMeasure>,
            g: &Graph,
            zeta: &Partition,
            eta: &Partition,
        ) -> f64;

        type ClusteringGenerator;
        pub fn NewClusteringGenerator() -> UniquePtr<ClusteringGenerator>;
        fn CMMakeContinuousBalancedClustering(
            gen: Pin<&mut ClusteringGenerator>,
            g: &Graph,
            k: u64,
        ) -> UniquePtr<Partition>;
        fn CMMakeNoncontinuousBalancedClustering(
            gen: Pin<&mut ClusteringGenerator>,
            g: &Graph,
            k: u64,
        ) -> UniquePtr<Partition>;
        fn CMMakeOneClustering(
            gen: Pin<&mut ClusteringGenerator>,
            g: &Graph,
        ) -> UniquePtr<Partition>;
        fn CMMakeRandomClustering(
            gen: Pin<&mut ClusteringGenerator>,
            g: &Graph,
            k: u64,
        ) -> UniquePtr<Partition>;
        fn CMMakeSingletonClustering(
            gen: Pin<&mut ClusteringGenerator>,
            g: &Graph,
        ) -> UniquePtr<Partition>;

        type CoverF1Similarity;
        fn getWeightedAverage(self: &CoverF1Similarity) -> f64;
        fn getUnweightedAverage(self: &CoverF1Similarity) -> f64;
        fn getMaximumValue(self: &CoverF1Similarity) -> f64;
        fn getMinimumValue(self: &CoverF1Similarity) -> f64;
        fn getValue(self: &CoverF1Similarity, i: u64) -> f64;
        fn CoverF1SimilarityGetValues(e: &CoverF1Similarity) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &CoverF1Similarity) -> bool;
        fn run(self: Pin<&mut CoverF1Similarity>) -> Result<()>;
        fn hasFinished(self: &CoverF1Similarity) -> bool;

        fn NewCoverF1Similarity(
            g: &Graph,
            c: &Cover,
            reference: &Cover,
        ) -> UniquePtr<CoverF1Similarity>;

        type CoverHubDominance;
        fn getWeightedAverage(self: &CoverHubDominance) -> f64;
        fn getUnweightedAverage(self: &CoverHubDominance) -> f64;
        fn getMaximumValue(self: &CoverHubDominance) -> f64;
        fn getMinimumValue(self: &CoverHubDominance) -> f64;
        fn getValue(self: &CoverHubDominance, i: u64) -> f64;
        fn CoverHubDominanceGetValues(e: &CoverHubDominance) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &CoverHubDominance) -> bool;
        fn run(self: Pin<&mut CoverHubDominance>) -> Result<()>;
        fn hasFinished(self: &CoverHubDominance) -> bool;

        fn NewCoverHubDominance(g: &Graph, c: &Cover) -> UniquePtr<CoverHubDominance>;

        type Coverage;
        fn NewCoverage() -> UniquePtr<Coverage>;
        fn getQuality(self: Pin<&mut Coverage>, p: &Partition, g: &Graph) -> f64;

        type CutClustering;
        fn NewCutClustering(g: &Graph, alpha: f64) -> UniquePtr<CutClustering>;
        // omitted CutClustering::getClusterHierarchy
        fn run(self: Pin<&mut CutClustering>) -> Result<()>;
        fn hasFinished(self: &CutClustering) -> bool;
        fn CutClusteringGetPartition(a: Pin<&mut CutClustering>) -> UniquePtr<Partition>;
    }
    #[namespace = "NetworKit::GraphTools"]
    unsafe extern "C++" {

        // GRAPH TOOLS

        fn append(g: Pin<&mut Graph>, g1: &Graph);
        fn augmentGraph(g: Pin<&mut Graph>) -> u64;
        fn GTCopyNodes(g: &Graph) -> UniquePtr<Graph>;
        fn GTCreateAugmentedGraph(g: &Graph, root: &mut u64) -> UniquePtr<Graph>;
        fn density(g: &Graph) -> f64;
        // getCompactedGraph, getContinuousNodeIds, getRandomContinuousNodeIds merged into one function
        fn GTGetCompactedGraph(g: &Graph, random: bool) -> UniquePtr<Graph>;
        fn GTVolume(g: &Graph, nodes: &[u64]) -> f64;
        fn GTInVolume(g: &Graph, nodes: &[u64]) -> f64;
        fn maxDegree(g: &Graph) -> u64;
        fn maxInDegree(g: &Graph) -> u64;
        fn maxWeightedDegree(g: &Graph) -> f64;
        fn maxWeightedInDegree(g: &Graph) -> f64;
        fn merge(g: Pin<&mut Graph>, g1: &Graph);
        fn GTRandomEdge(g: &Graph, uniform: bool, src: &mut u64, dst: &mut u64);
        fn GTRandomEdges(g: &Graph, n: u64, src: &mut Vec<u64>, dst: &mut Vec<u64>);
        unsafe fn randomNeighbor(g: &Graph, u: u64) -> u64;
        fn randomNode(g: &Graph) -> u64;
        fn GTRandomNodes(g: &Graph, n: u64) -> UniquePtr<CxxVector<u64>>;
        fn GTRemoveEdgesFromIsolatedSet(g: Pin<&mut Graph>, nodes: &[u64]);
        fn GTSize(g: &Graph, n_nodes: &mut u64, n_edges: &mut u64);
        fn sortEdgesByWeight(g: Pin<&mut Graph>, descending: bool);
        fn GTSubgraphAndNeighborsFromNodes(
            g: &Graph,
            nodes: &[u64],
            include_out_neighbors: bool,
            include_in_neighbors: bool,
        ) -> UniquePtr<Graph>;
        fn GTSubgraphFromNodes(g: &Graph, nodes: &[u64]) -> UniquePtr<Graph>;
        fn GTToUndirected(g: &Graph) -> UniquePtr<Graph>;
        fn GTToUnweighted(g: &Graph) -> UniquePtr<Graph>;
        fn GTToWeighted(g: &Graph) -> UniquePtr<Graph>;
        fn GTTopologicalSort(g: &Graph) -> UniquePtr<CxxVector<u64>>;
        fn GTTranspose(g: &Graph) -> Result<UniquePtr<Graph>>;
    }
}
