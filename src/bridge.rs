extern crate openmp_sys;

pub(crate) use ffi::*;

#[cxx::bridge]
mod ffi {

    #[namespace = "NetworKit"]
    unsafe extern "C++" {
        include!("bridge.h");

        // ---- GRAPH ----

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

        // ---- GRAPH BUILDER ----

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

        // ---- PARTITION ----

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

        // ---- COVER ----

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

        // ---- COMMUNITY ----

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
        // for CutClustering::getClusterHierarchy
        type HierarchyIter;
        fn isAtEnd(self: &HierarchyIter) -> bool;
        fn advance(self: Pin<&mut HierarchyIter>);
        fn curKey(self: &HierarchyIter) -> f64;
        fn curVal(self: &HierarchyIter) -> UniquePtr<Partition>;
        fn CutClusteringGetClusterHierarchy(g: &Graph) -> UniquePtr<HierarchyIter>;

        fn run(self: Pin<&mut CutClustering>) -> Result<()>;
        fn hasFinished(self: &CutClustering) -> bool;
        fn CutClusteringGetPartition(a: Pin<&mut CutClustering>) -> UniquePtr<Partition>;

        type EdgeCut;
        fn NewEdgeCut() -> UniquePtr<EdgeCut>;
        fn getQuality(self: Pin<&mut EdgeCut>, p: &Partition, g: &Graph) -> f64;

        type GraphStructuralRandMeasure;
        fn NewGraphStructuralRandMeasure() -> UniquePtr<GraphStructuralRandMeasure>;
        fn getDissimilarity(
            self: Pin<&mut GraphStructuralRandMeasure>,
            g: &Graph,
            first: &Partition,
            second: &Partition,
        ) -> f64;

        type HubDominance;
        fn NewHubDominance() -> UniquePtr<HubDominance>;
        fn getQuality(self: Pin<&mut HubDominance>, p: &Partition, g: &Graph) -> f64;
        #[rust_name = "getQualityForCover"]
        fn getQuality(self: Pin<&mut HubDominance>, p: &Cover, g: &Graph) -> f64;

        type IntrapartitionDensity;
        fn NewIntrapartitionDensity(g: &Graph, p: &Partition) -> UniquePtr<IntrapartitionDensity>;
        fn getGlobal(self: &IntrapartitionDensity) -> f64;
        fn getWeightedAverage(self: &IntrapartitionDensity) -> f64;
        fn getUnweightedAverage(self: &IntrapartitionDensity) -> f64;
        fn getMaximumValue(self: &IntrapartitionDensity) -> f64;
        fn getMinimumValue(self: &IntrapartitionDensity) -> f64;
        fn getValue(self: &IntrapartitionDensity, i: u64) -> f64;
        fn IntrapartitionDensityGetValues(e: &IntrapartitionDensity) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &IntrapartitionDensity) -> bool;
        fn run(self: Pin<&mut IntrapartitionDensity>) -> Result<()>;
        fn hasFinished(self: &IntrapartitionDensity) -> bool;

        type IsolatedInterpartitionConductance;
        fn NewIsolatedInterpartitionConductance(
            g: &Graph,
            p: &Partition,
        ) -> UniquePtr<IsolatedInterpartitionConductance>;
        fn getWeightedAverage(self: &IsolatedInterpartitionConductance) -> f64;
        fn getUnweightedAverage(self: &IsolatedInterpartitionConductance) -> f64;
        fn getMaximumValue(self: &IsolatedInterpartitionConductance) -> f64;
        fn getMinimumValue(self: &IsolatedInterpartitionConductance) -> f64;
        fn getValue(self: &IsolatedInterpartitionConductance, i: u64) -> f64;
        fn IsolatedInterpartitionConductanceGetValues(
            e: &IsolatedInterpartitionConductance,
        ) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &IsolatedInterpartitionConductance) -> bool;
        fn run(self: Pin<&mut IsolatedInterpartitionConductance>) -> Result<()>;
        fn hasFinished(self: &IsolatedInterpartitionConductance) -> bool;

        type IsolatedInterpartitionExpansion;
        fn NewIsolatedInterpartitionExpansion(
            g: &Graph,
            p: &Partition,
        ) -> UniquePtr<IsolatedInterpartitionExpansion>;
        fn getWeightedAverage(self: &IsolatedInterpartitionExpansion) -> f64;
        fn getUnweightedAverage(self: &IsolatedInterpartitionExpansion) -> f64;
        fn getMaximumValue(self: &IsolatedInterpartitionExpansion) -> f64;
        fn getMinimumValue(self: &IsolatedInterpartitionExpansion) -> f64;
        fn getValue(self: &IsolatedInterpartitionExpansion, i: u64) -> f64;
        fn IsolatedInterpartitionExpansionGetValues(
            e: &IsolatedInterpartitionExpansion,
        ) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &IsolatedInterpartitionExpansion) -> bool;
        fn run(self: Pin<&mut IsolatedInterpartitionExpansion>) -> Result<()>;
        fn hasFinished(self: &IsolatedInterpartitionExpansion) -> bool;

        type JaccardMeasure;
        pub fn NewJaccardMeasure() -> UniquePtr<JaccardMeasure>;
        pub fn getDissimilarity(
            self: Pin<&mut JaccardMeasure>,
            g: &Graph,
            zeta: &Partition,
            eta: &Partition,
        ) -> f64;

        type LFM;
        pub fn NewLFM(g: &Graph, scd: Pin<&mut SelectiveCommunityDetector>) -> UniquePtr<LFM>;
        fn LFMGetCover(algo: &LFM) -> UniquePtr<Cover>;
        fn run(self: Pin<&mut LFM>) -> Result<()>;
        fn hasFinished(self: &LFM) -> bool;

        type LPDegreeOrdered;
        fn NewLPDegreeOrdered(g: &Graph) -> UniquePtr<LPDegreeOrdered>;
        fn run(self: Pin<&mut LPDegreeOrdered>) -> Result<()>;
        fn hasFinished(self: &LPDegreeOrdered) -> bool;
        fn LPDegreeOrderedGetPartition(a: Pin<&mut LPDegreeOrdered>) -> UniquePtr<Partition>;
        fn numberOfIterations(self: Pin<&mut LPDegreeOrdered>) -> u64;

        type LouvainMapEquation;
        fn NewLouvainMapEquation(
            g: &Graph,
            hierarchical: bool,
            max_iterations: u64,
            parallelization_strategy: &str,
        ) -> UniquePtr<LouvainMapEquation>;
        fn run(self: Pin<&mut LouvainMapEquation>) -> Result<()>;
        fn hasFinished(self: &LouvainMapEquation) -> bool;
        fn LouvainMapEquationGetPartition(a: Pin<&mut LouvainMapEquation>) -> UniquePtr<Partition>;

        type Modularity;
        fn NewModularity() -> UniquePtr<Modularity>;
        fn getQuality(self: Pin<&mut Modularity>, p: &Partition, g: &Graph) -> f64;

        type NMIDistance;
        fn NewNMIDistance() -> UniquePtr<NMIDistance>;

        pub fn getDissimilarity(
            self: Pin<&mut NMIDistance>,
            g: &Graph,
            zeta: &Partition,
            eta: &Partition,
        ) -> f64;

        type NodeStructuralRandMeasure;
        fn NewNodeStructuralRandMeasure() -> UniquePtr<NodeStructuralRandMeasure>;

        pub fn getDissimilarity(
            self: Pin<&mut NodeStructuralRandMeasure>,
            g: &Graph,
            zeta: &Partition,
            eta: &Partition,
        ) -> f64;

        type OverlappingNMIDistance;
        pub fn NewOverlappingNMIDistance(normalization: u8) -> UniquePtr<OverlappingNMIDistance>;
        pub fn getDissimilarity(
            self: Pin<&mut OverlappingNMIDistance>,
            g: &Graph,
            zeta: &Partition,
            eta: &Partition,
        ) -> f64;
        #[rust_name = "getDissimilarityForCover"]
        pub fn getDissimilarity(
            self: Pin<&mut OverlappingNMIDistance>,
            g: &Graph,
            zeta: &Cover,
            eta: &Cover,
        ) -> f64;

        type PLM;
        pub fn NewPLM(
            g: &Graph,
            refine: bool,
            gamma: f64,
            par: &str,
            max_iter: u64,
            turbo: bool,
            recurse: bool,
        ) -> Result<UniquePtr<PLM>>;
        fn PLMCoarsen(g: &Graph, zeta: &Partition, mapping: &mut Vec<u64>) -> UniquePtr<Graph>;
        fn PLMProlong(
            g: &Graph,
            zeta_coarse: &Partition,
            g_fine: &Graph,
            node_to_meta_node: &[u64],
        ) -> UniquePtr<Partition>;
        fn PLMGetPartition(a: Pin<&mut PLM>) -> UniquePtr<Partition>;
        fn run(self: Pin<&mut PLM>) -> Result<()>;
        fn hasFinished(self: &PLM) -> bool;

        type PLP;
        fn NewPLP(g: &Graph, theta: u64, max_iterations: u64) -> UniquePtr<PLP>;
        fn run(self: Pin<&mut PLP>) -> Result<()>;
        fn hasFinished(self: &PLP) -> bool;
        fn numberOfIterations(self: Pin<&mut PLP>) -> u64;
        fn PLPGetPartition(a: Pin<&mut PLP>) -> UniquePtr<Partition>;

        type ParallelLeiden;
        fn NewParallelLeiden(
            g: &Graph,
            iterations: u64,
            randomize: bool,
            gamma: f64,
        ) -> UniquePtr<ParallelLeiden>;
        fn run(self: Pin<&mut ParallelLeiden>) -> Result<()>;
        fn hasFinished(self: &ParallelLeiden) -> bool;
        fn ParallelLeidenGetPartition(a: Pin<&mut ParallelLeiden>) -> UniquePtr<Partition>;

        type PartitionFragmentation;
        fn NewPartitionFragmentation(g: &Graph, p: &Partition)
            -> UniquePtr<PartitionFragmentation>;
        fn getWeightedAverage(self: &PartitionFragmentation) -> f64;
        fn getUnweightedAverage(self: &PartitionFragmentation) -> f64;
        fn getMaximumValue(self: &PartitionFragmentation) -> f64;
        fn getMinimumValue(self: &PartitionFragmentation) -> f64;
        fn getValue(self: &PartitionFragmentation, i: u64) -> f64;
        fn PartitionFragmentationGetValues(e: &PartitionFragmentation)
            -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &PartitionFragmentation) -> bool;
        fn run(self: Pin<&mut PartitionFragmentation>) -> Result<()>;
        fn hasFinished(self: &PartitionFragmentation) -> bool;

        type PartitionHubDominance;
        fn NewPartitionHubDominance(g: &Graph, p: &Partition) -> UniquePtr<PartitionHubDominance>;
        fn getWeightedAverage(self: &PartitionHubDominance) -> f64;
        fn getUnweightedAverage(self: &PartitionHubDominance) -> f64;
        fn getMaximumValue(self: &PartitionHubDominance) -> f64;
        fn getMinimumValue(self: &PartitionHubDominance) -> f64;
        fn getValue(self: &PartitionHubDominance, i: u64) -> f64;
        fn PartitionHubDominanceGetValues(e: &PartitionHubDominance) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &PartitionHubDominance) -> bool;
        fn run(self: Pin<&mut PartitionHubDominance>) -> Result<()>;
        fn hasFinished(self: &PartitionHubDominance) -> bool;

        type PartitionIntersection;
        fn NewPartitionIntersection() -> UniquePtr<PartitionIntersection>;
        fn PartitionIntersectionCalculate(
            algo: Pin<&mut PartitionIntersection>,
            zeta: &Partition,
            eta: &Partition,
        ) -> UniquePtr<Partition>;

        type StablePartitionNodes;
        fn NewStablePartitionNodes(g: &Graph, p: &Partition) -> UniquePtr<StablePartitionNodes>;
        fn getWeightedAverage(self: &StablePartitionNodes) -> f64;
        fn getUnweightedAverage(self: &StablePartitionNodes) -> f64;
        fn getMaximumValue(self: &StablePartitionNodes) -> f64;
        fn getMinimumValue(self: &StablePartitionNodes) -> f64;
        fn getValue(self: &StablePartitionNodes, i: u64) -> f64;
        fn StablePartitionNodesGetValues(e: &StablePartitionNodes) -> UniquePtr<CxxVector<f64>>;
        fn isSmallBetter(self: &StablePartitionNodes) -> bool;
        fn run(self: Pin<&mut StablePartitionNodes>) -> Result<()>;
        fn hasFinished(self: &StablePartitionNodes) -> bool;
        fn isStable(self: &StablePartitionNodes, u: u64) -> bool;

        // ---- SCD ----

        type ApproximatePageRank;
        fn NewApproximatePageRank(
            g: &Graph,
            alpha: f64,
            epsilon: f64,
        ) -> UniquePtr<ApproximatePageRank>;

        fn ApproximatePageRankRun(
            algo: Pin<&mut ApproximatePageRank>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );

        type SelectiveCommunityDetector;

        type CliqueDetect;
        fn NewCliqueDetect(g: &Graph) -> UniquePtr<CliqueDetect>;
        fn CliqueDetectRun(
            algo: Pin<&mut CliqueDetect>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn CliqueDetectExpandOneCommunity(
            algo: Pin<&mut CliqueDetect>,
            seeds: &[u64],
            ret: &mut Vec<u64>,
        );
        fn CliqueDetectAsBase(
            algo: UniquePtr<CliqueDetect>,
        ) -> UniquePtr<SelectiveCommunityDetector>;

        type CombinedSCD;
        fn NewCombinedSCD(
            g: &Graph,
            first: Pin<&mut SelectiveCommunityDetector>,
            second: Pin<&mut SelectiveCommunityDetector>,
        ) -> UniquePtr<CombinedSCD>;
        fn CombinedSCDRun(
            algo: Pin<&mut CombinedSCD>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn CombinedSCDExpandOneCommunity(
            algo: Pin<&mut CombinedSCD>,
            seeds: &[u64],
            ret: &mut Vec<u64>,
        );
        fn CombinedSCDAsBase(algo: UniquePtr<CombinedSCD>)
            -> UniquePtr<SelectiveCommunityDetector>;

        type GCE;
        fn NewGCE(g: &Graph, q: &str) -> UniquePtr<GCE>;
        fn GCERun(algo: Pin<&mut GCE>, seeds: &[u64], ks: &mut Vec<u64>, vs: &mut Vec<u64>);
        fn GCEExpandOneCommunity(algo: Pin<&mut GCE>, seeds: &[u64], ret: &mut Vec<u64>);
        fn GCEAsBase(algo: UniquePtr<GCE>) -> UniquePtr<SelectiveCommunityDetector>;

        type LFMLocal;
        fn NewLFMLocal(g: &Graph, alpha: f64) -> UniquePtr<LFMLocal>;
        fn LFMLocalRun(
            algo: Pin<&mut LFMLocal>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn LFMLocalExpandOneCommunity(algo: Pin<&mut LFMLocal>, seeds: &[u64], ret: &mut Vec<u64>);
        fn LFMLocalAsBase(algo: UniquePtr<LFMLocal>) -> UniquePtr<SelectiveCommunityDetector>;

        type LocalT;
        fn NewLocalT(g: &Graph) -> UniquePtr<LocalT>;
        fn LocalTRun(algo: Pin<&mut LocalT>, seeds: &[u64], ks: &mut Vec<u64>, vs: &mut Vec<u64>);
        fn LocalTExpandOneCommunity(algo: Pin<&mut LocalT>, seeds: &[u64], ret: &mut Vec<u64>);
        fn LocalTAsBase(algo: UniquePtr<LocalT>) -> UniquePtr<SelectiveCommunityDetector>;

        type LocalTightnessExpansion;
        fn NewLocalTightnessExpansion(g: &Graph, alpha: f64) -> UniquePtr<LocalTightnessExpansion>;
        fn LocalTightnessExpansionRun(
            algo: Pin<&mut LocalTightnessExpansion>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn LocalTightnessExpansionExpandOneCommunity(
            algo: Pin<&mut LocalTightnessExpansion>,
            seeds: &[u64],
            ret: &mut Vec<u64>,
        );
        fn LocalTightnessExpansionAsBase(
            algo: UniquePtr<LocalTightnessExpansion>,
        ) -> UniquePtr<SelectiveCommunityDetector>;

        type PageRankNibble;
        fn NewPageRankNibble(g: &Graph, alpha: f64, epsilon: f64) -> UniquePtr<PageRankNibble>;
        fn PageRankNibbleRun(
            algo: Pin<&mut PageRankNibble>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn PageRankNibbleExpandOneCommunity(
            algo: Pin<&mut PageRankNibble>,
            seeds: &[u64],
            ret: &mut Vec<u64>,
        );
        fn PageRankNibbleAsBase(
            algo: UniquePtr<PageRankNibble>,
        ) -> UniquePtr<SelectiveCommunityDetector>;

        type RandomBFS;
        fn NewRandomBFS(g: &Graph, c: &Cover) -> UniquePtr<RandomBFS>;
        fn RandomBFSRun(
            algo: Pin<&mut RandomBFS>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn RandomBFSExpandOneCommunity(
            algo: Pin<&mut RandomBFS>,
            seeds: &[u64],
            ret: &mut Vec<u64>,
        );
        fn RandomBFSAsBase(algo: UniquePtr<RandomBFS>) -> UniquePtr<SelectiveCommunityDetector>;

        type SCDGroundTruthComparison;
        fn NewSCDGroundTruthComparison(
            g: &Graph,
            ground_truth: &Cover,
            ks: &[u64],
            vs: &[u64],
            ignore_seeds: bool,
        ) -> UniquePtr<SCDGroundTruthComparison>;

        fn SCDGroundTruthComparisonGetIndividualJaccard(
            algo: &SCDGroundTruthComparison,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );

        fn SCDGroundTruthComparisonGetIndividualPrecision(
            algo: &SCDGroundTruthComparison,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );

        fn SCDGroundTruthComparisonGetIndividualRecall(
            algo: &SCDGroundTruthComparison,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );

        fn SCDGroundTruthComparisonGetIndividualF1(
            algo: &SCDGroundTruthComparison,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );

        fn getAverageJaccard(self: &SCDGroundTruthComparison) -> f64;
        fn getAverageF1(self: &SCDGroundTruthComparison) -> f64;
        fn getAveragePrecision(self: &SCDGroundTruthComparison) -> f64;
        fn getAverageRecall(self: &SCDGroundTruthComparison) -> f64;

        fn run(self: Pin<&mut SCDGroundTruthComparison>) -> Result<()>;
        fn hasFinished(self: &SCDGroundTruthComparison) -> bool;

        type SetConductance;
        fn NewSetConductance(g: &Graph, community: &[u64]) -> UniquePtr<SetConductance>;
        fn getConductance(self: &SetConductance) -> f64;
        fn run(self: Pin<&mut SetConductance>) -> Result<()>;
        fn hasFinished(self: &SetConductance) -> bool;

        type TCE;
        fn NewTCE(g: &Graph, refine: bool, use_jaccard: bool) -> UniquePtr<TCE>;
        fn TCERun(algo: Pin<&mut TCE>, seeds: &[u64], ks: &mut Vec<u64>, vs: &mut Vec<u64>);
        fn TCEExpandOneCommunity(algo: Pin<&mut TCE>, seeds: &[u64], ret: &mut Vec<u64>);
        fn TCEAsBase(algo: UniquePtr<TCE>) -> UniquePtr<SelectiveCommunityDetector>;

        type TwoPhaseL;
        fn NewTwoPhaseL(g: &Graph) -> UniquePtr<TwoPhaseL>;
        fn TwoPhaseLRun(
            algo: Pin<&mut TwoPhaseL>,
            seeds: &[u64],
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn TwoPhaseLExpandOneCommunity(
            algo: Pin<&mut TwoPhaseL>,
            seeds: &[u64],
            ret: &mut Vec<u64>,
        );
        fn TwoPhaseLAsBase(algo: UniquePtr<TwoPhaseL>) -> UniquePtr<SelectiveCommunityDetector>;

        // ---- COARSENING ----

        type ParallelPartitionCoarsening;
        fn NewParallelPartitionCoarsening(
            g: &Graph,
            zeta: &Partition,
            parallel: bool,
        ) -> UniquePtr<ParallelPartitionCoarsening>;

        fn run(self: Pin<&mut ParallelPartitionCoarsening>) -> Result<()>;
        fn hasFinished(self: &ParallelPartitionCoarsening) -> bool;
        fn ParallelPartitionCoarseningGetCoarseGraph(
            algo: &ParallelPartitionCoarsening,
        ) -> UniquePtr<Graph>;
        fn ParallelPartitionCoarseningGetFineToCoarseNodeMapping(
            algo: &ParallelPartitionCoarsening,
        ) -> UniquePtr<CxxVector<u64>>;

        // ---- CLIQUE ----

        type MaximalCliques;
        fn NewMaximalCliques(g: &Graph, maximum_only: bool) -> UniquePtr<MaximalCliques>;
        fn MaximalCliquesGetCliques(
            algo: Pin<&mut MaximalCliques>,
            cliques: &mut Vec<u64>,
            nodes: &mut Vec<u64>,
        );
        fn run(self: Pin<&mut MaximalCliques>) -> Result<()>;
        fn hasFinished(self: &MaximalCliques) -> bool;

        // ---- CENTRALITY ----

        type ApproxBetweenness;
        fn NewApproxBetweenness(
            g: &Graph,
            epsilon: f64,
            delta: f64,
            universal_constant: f64,
        ) -> UniquePtr<ApproxBetweenness>;
        fn run(self: Pin<&mut ApproxBetweenness>) -> Result<()>;
        fn hasFinished(self: &ApproxBetweenness) -> bool;
        fn centralization(self: Pin<&mut ApproxBetweenness>) -> f64;
        fn maximum(self: Pin<&mut ApproxBetweenness>) -> f64;
        fn score(self: Pin<&mut ApproxBetweenness>, node: u64) -> f64;
        fn ApproxBetweennessRanking(
            algo: Pin<&mut ApproxBetweenness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn ApproxBetweennessScores(algo: Pin<&mut ApproxBetweenness>) -> UniquePtr<CxxVector<f64>>;

        type ApproxCloseness;
        fn NewApproxCloseness(
            g: &Graph,
            n_samples: u64,
            epsilon: f64,
            normalized: bool,
            t: u8,
        ) -> UniquePtr<ApproxCloseness>;
        fn run(self: Pin<&mut ApproxCloseness>) -> Result<()>;
        fn hasFinished(self: &ApproxCloseness) -> bool;
        fn centralization(self: Pin<&mut ApproxCloseness>) -> f64;
        fn maximum(self: Pin<&mut ApproxCloseness>) -> f64;
        fn score(self: Pin<&mut ApproxCloseness>, node: u64) -> f64;
        fn ApproxClosenessRanking(
            algo: Pin<&mut ApproxCloseness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn ApproxClosenessScores(algo: Pin<&mut ApproxCloseness>) -> UniquePtr<CxxVector<f64>>;
        fn ApproxClosenessGetSquareErrorEstimates(
            algo: Pin<&mut ApproxCloseness>,
        ) -> UniquePtr<CxxVector<f64>>;

        type ApproxElectricalCloseness;
        fn NewApproxElectricalCloseness(
            g: &Graph,
            epsilon: f64,
            kappa: f64,
        ) -> UniquePtr<ApproxElectricalCloseness>;
        fn run(self: Pin<&mut ApproxElectricalCloseness>) -> Result<()>;
        fn hasFinished(self: &ApproxElectricalCloseness) -> bool;
        fn centralization(self: Pin<&mut ApproxElectricalCloseness>) -> f64;
        fn maximum(self: Pin<&mut ApproxElectricalCloseness>) -> f64;
        fn score(self: Pin<&mut ApproxElectricalCloseness>, node: u64) -> f64;
        fn ApproxElectricalClosenessRanking(
            algo: Pin<&mut ApproxElectricalCloseness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn ApproxElectricalClosenessScores(
            algo: Pin<&mut ApproxElectricalCloseness>,
        ) -> UniquePtr<CxxVector<f64>>;
        fn ApproxElectricalClosenessComputeExactDiagonal(
            algo: &ApproxElectricalCloseness,
            tol: f64,
        ) -> UniquePtr<CxxVector<f64>>;
        fn ApproxElectricalClosenessGetDiagonal(
            algo: &ApproxElectricalCloseness,
        ) -> UniquePtr<CxxVector<f64>>;

        type ApproxGroupBetweenness;
        fn NewApproxGroupBetweenness(
            g: &Graph,
            group_size: u64,
            epsilon: f64,
        ) -> UniquePtr<ApproxGroupBetweenness>;
        fn run(self: Pin<&mut ApproxGroupBetweenness>) -> Result<()>;
        fn hasFinished(self: &ApproxGroupBetweenness) -> bool;
        fn ApproxGroupBetweennessGroupMaxBetweenness(
            algo: &ApproxGroupBetweenness,
        ) -> UniquePtr<CxxVector<u64>>;
        fn ApproxGroupBetweennessScoreOfGroup(
            algo: &ApproxGroupBetweenness,
            nodes: &[u64],
            normalized: bool,
        ) -> UniquePtr<CxxVector<u64>>;

        type ApproxSpanningEdge;
        fn NewApproxSpanningEdge(g: &Graph, epsilon: f64) -> UniquePtr<ApproxSpanningEdge>;
        fn ApproxSpanningEdgeScores(algo: &ApproxSpanningEdge) -> UniquePtr<CxxVector<f64>>;
        fn run(self: Pin<&mut ApproxSpanningEdge>) -> Result<()>;
        fn hasFinished(self: &ApproxSpanningEdge) -> bool;

        type Betweenness;
        fn NewBetweenness(
            g: &Graph,
            normalized: bool,
            compute_edge_centrality: bool,
        ) -> UniquePtr<Betweenness>;
        fn run(self: Pin<&mut Betweenness>) -> Result<()>;
        fn hasFinished(self: &Betweenness) -> bool;
        fn centralization(self: Pin<&mut Betweenness>) -> f64;
        fn maximum(self: Pin<&mut Betweenness>) -> f64;
        fn score(self: Pin<&mut Betweenness>, node: u64) -> f64;
        fn BetweennessRanking(algo: Pin<&mut Betweenness>, ks: &mut Vec<u64>, vs: &mut Vec<f64>);
        fn BetweennessScores(algo: Pin<&mut Betweenness>) -> UniquePtr<CxxVector<f64>>;
        fn BetweennessEdgeScores(algo: Pin<&mut Betweenness>) -> UniquePtr<CxxVector<f64>>;

        type Closeness;
        fn NewCloseness(g: &Graph, normalized: bool, variant: u8) -> UniquePtr<Closeness>;
        fn run(self: Pin<&mut Closeness>) -> Result<()>;
        fn hasFinished(self: &Closeness) -> bool;
        fn centralization(self: Pin<&mut Closeness>) -> f64;
        fn maximum(self: Pin<&mut Closeness>) -> f64;
        fn score(self: Pin<&mut Closeness>, node: u64) -> f64;
        fn ClosenessRanking(algo: Pin<&mut Closeness>, ks: &mut Vec<u64>, vs: &mut Vec<f64>);
        fn ClosenessScores(algo: Pin<&mut Closeness>) -> UniquePtr<CxxVector<f64>>;

        type CoreDecomposition;
        fn NewCoreDecomposition(
            g: &Graph,
            normalized: bool,
            enforce_bucket_queue_algorithm: bool,
            store_node_order: bool,
        ) -> UniquePtr<CoreDecomposition>;
        fn run(self: Pin<&mut CoreDecomposition>) -> Result<()>;
        fn hasFinished(self: &CoreDecomposition) -> bool;
        fn centralization(self: Pin<&mut CoreDecomposition>) -> f64;
        fn maximum(self: Pin<&mut CoreDecomposition>) -> f64;
        fn score(self: Pin<&mut CoreDecomposition>, node: u64) -> f64;
        fn CoreDecompositionRanking(
            algo: Pin<&mut CoreDecomposition>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn CoreDecompositionScores(algo: Pin<&mut CoreDecomposition>) -> UniquePtr<CxxVector<f64>>;
        fn maxCoreNumber(self: &CoreDecomposition) -> u64;
        fn CoreDecompositionGetCover(algo: &CoreDecomposition) -> UniquePtr<Cover>;
        fn CoreDecompositionGetPartition(algo: &CoreDecomposition) -> UniquePtr<Partition>;
        fn CoreDecompositionGetNodeOrder(algo: &CoreDecomposition) -> UniquePtr<CxxVector<u64>>;

        type DegreeCentrality;
        fn NewDegreeCentrality(
            g: &Graph,
            normalized: bool,
            out_deg: bool,
            ignore_self_loops: bool,
        ) -> UniquePtr<DegreeCentrality>;
        fn run(self: Pin<&mut DegreeCentrality>) -> Result<()>;
        fn hasFinished(self: &DegreeCentrality) -> bool;
        fn centralization(self: Pin<&mut DegreeCentrality>) -> f64;
        fn maximum(self: Pin<&mut DegreeCentrality>) -> f64;
        fn score(self: Pin<&mut DegreeCentrality>, node: u64) -> f64;
        fn DegreeCentralityRanking(
            algo: Pin<&mut DegreeCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn DegreeCentralityScores(algo: Pin<&mut DegreeCentrality>) -> UniquePtr<CxxVector<f64>>;

        type DynApproxBetweenness;
        fn NewDynApproxBetweenness(
            g: &Graph,
            epsilon: f64,
            delta: f64,
            store_predecessors: bool,
            universal_constant: f64,
        ) -> UniquePtr<DynApproxBetweenness>;
        fn run(self: Pin<&mut DynApproxBetweenness>) -> Result<()>;
        fn hasFinished(self: &DynApproxBetweenness) -> bool;
        fn centralization(self: Pin<&mut DynApproxBetweenness>) -> f64;
        fn maximum(self: Pin<&mut DynApproxBetweenness>) -> f64;
        fn score(self: Pin<&mut DynApproxBetweenness>, node: u64) -> f64;
        fn DynApproxBetweennessRanking(
            algo: Pin<&mut DynApproxBetweenness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn DynApproxBetweennessScores(
            algo: Pin<&mut DynApproxBetweenness>,
        ) -> UniquePtr<CxxVector<f64>>;
        fn getNumberOfSamples(self: &DynApproxBetweenness) -> u64;
        fn DynApproxBetweennessUpdate(
            algo: Pin<&mut DynApproxBetweenness>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        );
        fn DynApproxBetweennessUpdateBatch(
            algo: Pin<&mut DynApproxBetweenness>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );
    }
    #[namespace = "NetworKit::GraphTools"]
    unsafe extern "C++" {

        // ---- GRAPH TOOLS ----

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
    #[namespace = "NetworKit::GraphClusteringTools"]
    unsafe extern "C++" {
        fn MakeCommunicationGraph(g: &Graph, zeta: Pin<&mut Partition>) -> UniquePtr<Graph>;
        fn equalClusterings(zeta: &Partition, eta: &Partition, g: Pin<&mut Graph>) -> bool;
        fn getImbalance(zeta: &Partition) -> f32;
        fn isOneClustering(g: &Graph, zeta: &Partition) -> bool;
        fn isProperClustering(g: &Graph, zeta: &Partition) -> bool;
        fn isSingletonClustering(g: &Graph, zeta: &Partition) -> bool;
        fn weightedDegreeWithCluster(g: &Graph, zeta: &Partition, u: u64, cid: u64) -> u64;
    }
}
