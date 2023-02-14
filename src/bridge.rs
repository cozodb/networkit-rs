extern crate openmp_sys;

pub(crate) use ffi::*;

#[cxx::bridge]
mod ffi {

    #[namespace = "NetworKit"]
    unsafe extern "C++" {
        include!("bridge.h");

        fn MakeWeightVector(wt: &[f64]) -> UniquePtr<CxxVector<f64>>;
        fn MakeCountVector(ct: &[u64]) -> UniquePtr<CxxVector<u64>>;

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

        type DynBetweenness;
        fn NewDynBetweenness(g: &Graph) -> UniquePtr<DynBetweenness>;
        fn run(self: Pin<&mut DynBetweenness>) -> Result<()>;
        fn hasFinished(self: &DynBetweenness) -> bool;
        fn centralization(self: Pin<&mut DynBetweenness>) -> f64;
        fn maximum(self: Pin<&mut DynBetweenness>) -> f64;
        fn score(self: Pin<&mut DynBetweenness>, node: u64) -> f64;
        fn DynBetweennessRanking(
            algo: Pin<&mut DynBetweenness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn DynBetweennessScores(algo: Pin<&mut DynBetweenness>) -> UniquePtr<CxxVector<f64>>;
        fn DynBetweennessUpdate(algo: Pin<&mut DynBetweenness>, kind: u8, u: u64, v: u64, ew: f64);
        fn DynBetweennessUpdateBatch(
            algo: Pin<&mut DynBetweenness>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );

        type DynBetweennessOneNode;
        fn NewDynBetweennessOneNode(g: Pin<&mut Graph>, x: u64)
            -> UniquePtr<DynBetweennessOneNode>;
        fn DynBetweennessOneNodeUpdate(
            algo: Pin<&mut DynBetweennessOneNode>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        );
        fn DynBetweennessOneNodeUpdateBatch(
            algo: Pin<&mut DynBetweennessOneNode>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );
        fn DynBetweennessOneNodeComputeScore(
            algo: Pin<&mut DynBetweennessOneNode>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        ) -> f64;
        fn getDistance(self: Pin<&mut DynBetweennessOneNode>, u: u64, v: u64) -> f64;
        fn getSigma(self: Pin<&mut DynBetweennessOneNode>, u: u64, v: u64) -> f64;
        fn getSigmax(self: Pin<&mut DynBetweennessOneNode>, u: u64, v: u64) -> f64;
        fn getbcx(self: Pin<&mut DynBetweennessOneNode>) -> f64;
        fn run(self: Pin<&mut DynBetweennessOneNode>);

        type DynKatzCentrality;
        fn NewDynKatzCentrality(
            g: &Graph,
            k: u64,
            group_only: bool,
            tolerance: f64,
        ) -> UniquePtr<DynKatzCentrality>;
        fn areDistinguished(self: Pin<&mut DynKatzCentrality>, u: u64, v: u64) -> bool;
        fn bound(self: Pin<&mut DynKatzCentrality>, v: u64) -> f64;
        fn run(self: Pin<&mut DynKatzCentrality>) -> Result<()>;
        fn hasFinished(self: &DynKatzCentrality) -> bool;
        fn centralization(self: Pin<&mut DynKatzCentrality>) -> f64;
        fn maximum(self: Pin<&mut DynKatzCentrality>) -> f64;
        fn score(self: Pin<&mut DynKatzCentrality>, node: u64) -> f64;
        fn DynKatzCentralityRanking(
            algo: Pin<&mut DynKatzCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn DynKatzCentralityScores(algo: Pin<&mut DynKatzCentrality>) -> UniquePtr<CxxVector<f64>>;
        fn DynKatzCentralityUpdate(
            algo: Pin<&mut DynKatzCentrality>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        );
        fn DynKatzCentralityUpdateBatch(
            algo: Pin<&mut DynKatzCentrality>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );
        fn DynKatzCentralityTop(
            algo: Pin<&mut DynKatzCentrality>,
            n: u64,
        ) -> UniquePtr<CxxVector<u64>>;

        type DynTopHarmonicCloseness;
        fn NewDynTopHarmonicCloseness(
            g: &Graph,
            k: u64,
            use_bfs_bound: bool,
        ) -> UniquePtr<DynTopHarmonicCloseness>;
        fn run(self: Pin<&mut DynTopHarmonicCloseness>) -> Result<()>;
        fn hasFinished(self: &DynTopHarmonicCloseness) -> bool;
        fn DynTopHarmonicClosenessRanking(
            algo: Pin<&mut DynTopHarmonicCloseness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn DynTopHarmonicClosenessUpdate(
            algo: Pin<&mut DynTopHarmonicCloseness>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        );
        fn DynTopHarmonicClosenessUpdateBatch(
            algo: Pin<&mut DynTopHarmonicCloseness>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );
        fn reset(self: Pin<&mut DynTopHarmonicCloseness>);
        fn DynTopHarmonicClosenessTopkNodesList(
            algo: Pin<&mut DynTopHarmonicCloseness>,
            include_trail: bool,
        ) -> UniquePtr<CxxVector<u64>>;
        fn DynTopHarmonicClosenessTopkScoresList(
            algo: Pin<&mut DynTopHarmonicCloseness>,
            include_trail: bool,
        ) -> UniquePtr<CxxVector<f64>>;

        type EigenvectorCentrality;
        fn NewEigenvectorCentrality(g: &Graph, tol: f64) -> UniquePtr<EigenvectorCentrality>;
        fn run(self: Pin<&mut EigenvectorCentrality>) -> Result<()>;
        fn hasFinished(self: &EigenvectorCentrality) -> bool;
        fn centralization(self: Pin<&mut EigenvectorCentrality>) -> f64;
        fn maximum(self: Pin<&mut EigenvectorCentrality>) -> f64;
        fn score(self: Pin<&mut EigenvectorCentrality>, node: u64) -> f64;
        fn EigenvectorCentralityRanking(
            algo: Pin<&mut EigenvectorCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn EigenvectorCentralityScores(
            algo: Pin<&mut EigenvectorCentrality>,
        ) -> UniquePtr<CxxVector<f64>>;

        type EstimateBetweenness;
        fn NewEstimateBetweenness(
            g: &Graph,
            n_samples: u64,
            normalized: bool,
            parallel: bool,
        ) -> UniquePtr<EstimateBetweenness>;
        fn run(self: Pin<&mut EstimateBetweenness>) -> Result<()>;
        fn hasFinished(self: &EstimateBetweenness) -> bool;
        fn centralization(self: Pin<&mut EstimateBetweenness>) -> f64;
        fn maximum(self: Pin<&mut EstimateBetweenness>) -> f64;
        fn score(self: Pin<&mut EstimateBetweenness>, node: u64) -> f64;
        fn EstimateBetweennessRanking(
            algo: Pin<&mut EstimateBetweenness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn EstimateBetweennessScores(
            algo: Pin<&mut EstimateBetweenness>,
        ) -> UniquePtr<CxxVector<f64>>;

        type ForestCentrality;
        fn NewForestCentrality(
            g: &Graph,
            root: u64,
            epsilon: f64,
            kappa: f64,
        ) -> UniquePtr<ForestCentrality>;
        fn run(self: Pin<&mut ForestCentrality>) -> Result<()>;
        fn hasFinished(self: &ForestCentrality) -> bool;
        fn centralization(self: Pin<&mut ForestCentrality>) -> f64;
        fn maximum(self: Pin<&mut ForestCentrality>) -> f64;
        fn score(self: Pin<&mut ForestCentrality>, node: u64) -> f64;
        fn ForestCentralityRanking(
            algo: Pin<&mut ForestCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn ForestCentralityScores(algo: Pin<&mut ForestCentrality>) -> UniquePtr<CxxVector<f64>>;
        fn ForestCentralityGetDiagonal(algo: &ForestCentrality) -> UniquePtr<CxxVector<f64>>;
        fn getNumberOfSamples(self: &ForestCentrality) -> u64;

        type GedWalk;
        fn NewGedWalk(
            g: &Graph,
            k: u64,
            init_epsilon: f64,
            alpha: f64,
            bs: u8,
            gs: u8,
            spectral_delta: f64,
        ) -> UniquePtr<GedWalk>;
        fn run(self: Pin<&mut GedWalk>) -> Result<()>;
        fn hasFinished(self: &GedWalk) -> bool;
        fn getApproximateScore(self: &GedWalk) -> f64;
        fn GedWalkGroupMaxGedWalk(algo: &GedWalk) -> UniquePtr<CxxVector<u64>>;
        fn GedWalkScoreOfGroup(algo: Pin<&mut GedWalk>, group: &[u64], epsilon: f64) -> f64;

        type GroupCloseness;
        fn NewGroupCloseness(g: &Graph, k: u64, h: u64) -> UniquePtr<GroupCloseness>;
        fn run(self: Pin<&mut GroupCloseness>) -> Result<()>;
        fn hasFinished(self: &GroupCloseness) -> bool;
        fn GroupClosenessScoreOfGroup(algo: Pin<&mut GroupCloseness>, group: &[u64]) -> f64;
        fn GroupClosenessGroupMaxCloseness(
            algo: Pin<&mut GroupCloseness>,
        ) -> UniquePtr<CxxVector<u64>>;
        fn GroupClosenessComputeFarness(algo: &GroupCloseness, group: &[u64], h: u64) -> f64;

        type GroupClosenessGrowShrink;
        fn NewGroupClosenessGrowShrink(
            g: &Graph,
            group: &[u64],
            extended: bool,
            insertions: u64,
            max_iterations: u64,
        ) -> UniquePtr<GroupClosenessGrowShrink>;
        fn run(self: Pin<&mut GroupClosenessGrowShrink>) -> Result<()>;
        fn hasFinished(self: &GroupClosenessGrowShrink) -> bool;
        fn GroupClosenessGrowShrinkGroupMaxCloseness(
            algo: &GroupClosenessGrowShrink,
        ) -> UniquePtr<CxxVector<u64>>;
        fn numberOfIterations(self: &GroupClosenessGrowShrink) -> u64;

        type GroupClosenessLocalSearch;
        fn NewGroupClosenessLocalSearch(
            g: &Graph,
            group: &[u64],
            run_grow_shrink: bool,
            max_iterations: u64,
        ) -> UniquePtr<GroupClosenessLocalSearch>;
        fn run(self: Pin<&mut GroupClosenessLocalSearch>) -> Result<()>;
        fn hasFinished(self: &GroupClosenessLocalSearch) -> bool;
        fn GroupClosenessLocalSearchGroupMaxCloseness(
            algo: &GroupClosenessLocalSearch,
        ) -> UniquePtr<CxxVector<u64>>;
        fn numberOfIterations(self: &GroupClosenessLocalSearch) -> u64;

        type GroupClosenessLocalSwaps;
        fn NewGroupClosenessLocalSwaps(
            g: &Graph,
            group: &[u64],
            max_swaps: u64,
        ) -> UniquePtr<GroupClosenessLocalSwaps>;
        fn run(self: Pin<&mut GroupClosenessLocalSwaps>) -> Result<()>;
        fn hasFinished(self: &GroupClosenessLocalSwaps) -> bool;
        fn GroupClosenessLocalSwapsGroupMaxCloseness(
            algo: &GroupClosenessLocalSwaps,
        ) -> UniquePtr<CxxVector<u64>>;
        fn numberOfSwaps(self: &GroupClosenessLocalSwaps) -> u64;

        type GroupDegree;
        fn NewGroupDegree(g: &Graph, k: u64, count_group_nodes: bool) -> UniquePtr<GroupDegree>;
        fn run(self: Pin<&mut GroupDegree>) -> Result<()>;
        fn hasFinished(self: &GroupDegree) -> bool;
        fn GroupDegreeGroupMaxDegree(algo: Pin<&mut GroupDegree>) -> UniquePtr<CxxVector<u64>>;
        fn GroupDegreeScoreOfGroup(algo: &GroupDegree, group: &[u64]) -> f64;
        fn getScore(self: Pin<&mut GroupDegree>) -> u64;

        type GroupHarmonicCloseness;
        fn NewGroupHarmonicCloseness(g: &Graph, k: u64) -> UniquePtr<GroupHarmonicCloseness>;
        fn run(self: Pin<&mut GroupHarmonicCloseness>) -> Result<()>;
        fn hasFinished(self: &GroupHarmonicCloseness) -> bool;
        fn GroupHarmonicClosenessGroupMaxHarmonicCloseness(
            algo: Pin<&mut GroupHarmonicCloseness>,
        ) -> UniquePtr<CxxVector<u64>>;
        fn GroupHarmonicClosenessScoreOfGroup(g: &Graph, group: &[u64]) -> f64;

        type HarmonicCloseness;
        fn NewHarmonicCloseness(g: &Graph, normalized: bool) -> UniquePtr<HarmonicCloseness>;
        fn run(self: Pin<&mut HarmonicCloseness>) -> Result<()>;
        fn hasFinished(self: &HarmonicCloseness) -> bool;
        fn centralization(self: Pin<&mut HarmonicCloseness>) -> f64;
        fn maximum(self: Pin<&mut HarmonicCloseness>) -> f64;
        fn score(self: Pin<&mut HarmonicCloseness>, node: u64) -> f64;
        fn HarmonicClosenessRanking(
            algo: Pin<&mut HarmonicCloseness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn HarmonicClosenessScores(algo: Pin<&mut HarmonicCloseness>) -> UniquePtr<CxxVector<f64>>;

        type KPathCentrality;
        fn NewKPathCentrality(g: &Graph, alpha: f64, k: u64) -> UniquePtr<KPathCentrality>;
        fn run(self: Pin<&mut KPathCentrality>) -> Result<()>;
        fn hasFinished(self: &KPathCentrality) -> bool;
        fn centralization(self: Pin<&mut KPathCentrality>) -> f64;
        fn maximum(self: Pin<&mut KPathCentrality>) -> f64;
        fn score(self: Pin<&mut KPathCentrality>, node: u64) -> f64;
        fn KPathCentralityRanking(
            algo: Pin<&mut KPathCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn KPathCentralityScores(algo: Pin<&mut KPathCentrality>) -> UniquePtr<CxxVector<f64>>;

        type KadabraBetweenness;
        fn NewKadabraBetweenness(
            g: &Graph,
            err: f64,
            delta: f64,
            deterministic: bool,
            k: u64,
            union_sample: u64,
            start_factor: u64,
        ) -> UniquePtr<KadabraBetweenness>;
        fn run(self: Pin<&mut KadabraBetweenness>) -> Result<()>;
        fn hasFinished(self: &KadabraBetweenness) -> bool;
        fn KadabraBetweennessRanking(
            algo: Pin<&mut KadabraBetweenness>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn KadabraBetweennessScores(
            algo: Pin<&mut KadabraBetweenness>,
        ) -> UniquePtr<CxxVector<f64>>;
        fn getNumberOfIterations(self: &KadabraBetweenness) -> u64;
        fn getOmega(self: &KadabraBetweenness) -> f64;

        fn KadabraBetweennessTopkNodesList(
            algo: Pin<&mut KadabraBetweenness>,
        ) -> UniquePtr<CxxVector<u64>>;
        fn KadabraBetweennessTopkScoresList(
            algo: Pin<&mut KadabraBetweenness>,
        ) -> UniquePtr<CxxVector<f64>>;

        type KatzCentrality;
        fn NewKatzCentrality(
            g: &Graph,
            alpha: f64,
            beta: f64,
            tol: f64,
        ) -> UniquePtr<KatzCentrality>;
        fn run(self: Pin<&mut KatzCentrality>) -> Result<()>;
        fn hasFinished(self: &KatzCentrality) -> bool;
        fn centralization(self: Pin<&mut KatzCentrality>) -> f64;
        fn maximum(self: Pin<&mut KatzCentrality>) -> f64;
        fn score(self: Pin<&mut KatzCentrality>, node: u64) -> f64;
        fn KatzCentralityRanking(
            algo: Pin<&mut KatzCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn KatzCentralityScores(algo: Pin<&mut KatzCentrality>) -> UniquePtr<CxxVector<f64>>;
        fn KatzCentralitySetEdgeDirection(algo: Pin<&mut KatzCentrality>, is_out: bool);

        type LaplacianCentrality;
        fn NewLaplacianCentrality(g: &Graph, normalized: bool) -> UniquePtr<LaplacianCentrality>;
        fn run(self: Pin<&mut LaplacianCentrality>) -> Result<()>;
        fn hasFinished(self: &LaplacianCentrality) -> bool;
        fn centralization(self: Pin<&mut LaplacianCentrality>) -> f64;
        fn maximum(self: Pin<&mut LaplacianCentrality>) -> f64;
        fn score(self: Pin<&mut LaplacianCentrality>, node: u64) -> f64;
        fn LaplacianCentralityRanking(
            algo: Pin<&mut LaplacianCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn LaplacianCentralityScores(
            algo: Pin<&mut LaplacianCentrality>,
        ) -> UniquePtr<CxxVector<f64>>;

        type LocalClusteringCoefficient;
        fn NewLocalClusteringCoefficient(
            g: &Graph,
            turbo: bool,
        ) -> UniquePtr<LocalClusteringCoefficient>;
        fn run(self: Pin<&mut LocalClusteringCoefficient>) -> Result<()>;
        fn hasFinished(self: &LocalClusteringCoefficient) -> bool;
        fn centralization(self: Pin<&mut LocalClusteringCoefficient>) -> f64;
        fn maximum(self: Pin<&mut LocalClusteringCoefficient>) -> f64;
        fn score(self: Pin<&mut LocalClusteringCoefficient>, node: u64) -> f64;
        fn LocalClusteringCoefficientRanking(
            algo: Pin<&mut LocalClusteringCoefficient>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn LocalClusteringCoefficientScores(
            algo: Pin<&mut LocalClusteringCoefficient>,
        ) -> UniquePtr<CxxVector<f64>>;

        type LocalPartitionCoverage;
        fn NewLocalPartitionCoverage(
            g: &Graph,
            partition: &Partition,
        ) -> UniquePtr<LocalPartitionCoverage>;
        fn run(self: Pin<&mut LocalPartitionCoverage>) -> Result<()>;
        fn hasFinished(self: &LocalPartitionCoverage) -> bool;
        fn centralization(self: Pin<&mut LocalPartitionCoverage>) -> f64;
        fn maximum(self: Pin<&mut LocalPartitionCoverage>) -> f64;
        fn score(self: Pin<&mut LocalPartitionCoverage>, node: u64) -> f64;
        fn LocalPartitionCoverageRanking(
            algo: Pin<&mut LocalPartitionCoverage>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn LocalPartitionCoverageScores(
            algo: Pin<&mut LocalPartitionCoverage>,
        ) -> UniquePtr<CxxVector<f64>>;

        type LocalSquareClusteringCoefficient;
        fn NewLocalSquareClusteringCoefficient(
            g: &Graph,
        ) -> UniquePtr<LocalSquareClusteringCoefficient>;
        fn run(self: Pin<&mut LocalSquareClusteringCoefficient>) -> Result<()>;
        fn hasFinished(self: &LocalSquareClusteringCoefficient) -> bool;
        fn centralization(self: Pin<&mut LocalSquareClusteringCoefficient>) -> f64;
        fn maximum(self: Pin<&mut LocalSquareClusteringCoefficient>) -> f64;
        fn score(self: Pin<&mut LocalSquareClusteringCoefficient>, node: u64) -> f64;
        fn LocalSquareClusteringCoefficientRanking(
            algo: Pin<&mut LocalSquareClusteringCoefficient>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn LocalSquareClusteringCoefficientScores(
            algo: Pin<&mut LocalSquareClusteringCoefficient>,
        ) -> UniquePtr<CxxVector<f64>>;

        type PageRank;
        fn NewPageRank(
            g: &Graph,
            damp: f64,
            tol: f64,
            normalized: bool,
            distribute_sinks: bool,
        ) -> UniquePtr<PageRank>;
        fn run(self: Pin<&mut PageRank>) -> Result<()>;
        fn hasFinished(self: &PageRank) -> bool;
        fn centralization(self: Pin<&mut PageRank>) -> f64;
        fn maximum(self: Pin<&mut PageRank>) -> f64;
        fn score(self: Pin<&mut PageRank>, node: u64) -> f64;
        fn PageRankRanking(algo: Pin<&mut PageRank>, ks: &mut Vec<u64>, vs: &mut Vec<f64>);
        fn PageRankScores(algo: Pin<&mut PageRank>) -> UniquePtr<CxxVector<f64>>;
        fn numberOfIterations(self: &PageRank) -> u64;
        fn PageRankSetMaxIterations(algo: Pin<&mut PageRank>, max_iter: u64);
        fn PageRankSetNorm(algo: Pin<&mut PageRank>, norm: u8);

        type PermanenceCentrality;
        fn NewPermanenceCentrality(g: &Graph, p: &Partition) -> UniquePtr<PermanenceCentrality>;
        fn run(self: Pin<&mut PermanenceCentrality>) -> Result<()>;
        fn hasFinished(self: &PermanenceCentrality) -> bool;
        fn getPermanence(self: Pin<&mut PermanenceCentrality>, u: u64) -> f64;
        fn getIntraClustering(self: Pin<&mut PermanenceCentrality>, u: u64) -> f64;

        type Sfigality;
        fn NewSfigality(g: &Graph) -> UniquePtr<Sfigality>;
        fn run(self: Pin<&mut Sfigality>) -> Result<()>;
        fn hasFinished(self: &Sfigality) -> bool;
        fn centralization(self: Pin<&mut Sfigality>) -> f64;
        fn maximum(self: Pin<&mut Sfigality>) -> f64;
        fn score(self: Pin<&mut Sfigality>, node: u64) -> f64;
        fn SfigalityRanking(algo: Pin<&mut Sfigality>, ks: &mut Vec<u64>, vs: &mut Vec<f64>);
        fn SfigalityScores(algo: Pin<&mut Sfigality>) -> UniquePtr<CxxVector<f64>>;

        type SpanningEdgeCentrality;
        fn NewSpanningEdgeCentrality(g: &Graph, tol: f64) -> UniquePtr<SpanningEdgeCentrality>;
        fn run(self: Pin<&mut SpanningEdgeCentrality>) -> Result<()>;
        fn hasFinished(self: &SpanningEdgeCentrality) -> bool;
        fn centralization(self: Pin<&mut SpanningEdgeCentrality>) -> f64;
        fn maximum(self: Pin<&mut SpanningEdgeCentrality>) -> f64;
        fn score(self: Pin<&mut SpanningEdgeCentrality>, node: u64) -> f64;
        fn SpanningEdgeCentralityRanking(
            algo: Pin<&mut SpanningEdgeCentrality>,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );
        fn SpanningEdgeCentralityScores(
            algo: Pin<&mut SpanningEdgeCentrality>,
        ) -> UniquePtr<CxxVector<f64>>;
        fn runApproximation(self: Pin<&mut SpanningEdgeCentrality>);
        fn runParallelApproximation(self: Pin<&mut SpanningEdgeCentrality>);

        type TopCloseness;
        fn NewTopCloseness(
            g: &Graph,
            k: u64,
            first_heu: bool,
            sec_heu: bool,
        ) -> UniquePtr<TopCloseness>;
        fn run(self: Pin<&mut TopCloseness>) -> Result<()>;
        fn hasFinished(self: &TopCloseness) -> bool;
        fn TopClosenessTopkNodesList(
            algo: Pin<&mut TopCloseness>,
            include_trail: bool,
        ) -> UniquePtr<CxxVector<u64>>;
        fn TopClosenessTopkScoresList(
            algo: Pin<&mut TopCloseness>,
            include_trail: bool,
        ) -> UniquePtr<CxxVector<f64>>;
        fn TopClosenessRestrictTopKComputationToNodes(algo: Pin<&mut TopCloseness>, nodes: &[u64]);

        type TopHarmonicCloseness;
        fn NewTopHarmonicCloseness(
            g: &Graph,
            k: u64,
            use_nb_bound: bool,
        ) -> UniquePtr<TopHarmonicCloseness>;
        fn run(self: Pin<&mut TopHarmonicCloseness>) -> Result<()>;
        fn hasFinished(self: &TopHarmonicCloseness) -> bool;
        fn TopHarmonicClosenessTopkNodesList(
            algo: Pin<&mut TopHarmonicCloseness>,
            include_trail: bool,
        ) -> UniquePtr<CxxVector<u64>>;
        fn TopHarmonicClosenessTopkScoresList(
            algo: Pin<&mut TopHarmonicCloseness>,
            include_trail: bool,
        ) -> UniquePtr<CxxVector<f64>>;
        fn TopHarmonicClosenessRestrictTopKComputationToNodes(
            algo: Pin<&mut TopHarmonicCloseness>,
            nodes: &[u64],
        );

        // ---- COMPONENTS ----

        type BiconnectedComponents;
        fn NewBiconnectedComponents(g: &Graph) -> UniquePtr<BiconnectedComponents>;
        fn run(self: Pin<&mut BiconnectedComponents>) -> Result<()>;
        fn hasFinished(self: &BiconnectedComponents) -> bool;
        fn numberOfComponents(self: &BiconnectedComponents) -> u64;
        fn BiconnectedComponentsGetComponentSizes(
            algo: &BiconnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn BiconnectedComponentsGetComponents(
            algo: &BiconnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn BiconnectedComponentsGetComponentOfNode(
            algo: &BiconnectedComponents,
            u: u64,
            vs: &mut Vec<u64>,
        );

        type ConnectedComponents;
        fn NewConnectedComponents(g: &Graph) -> UniquePtr<ConnectedComponents>;
        fn run(self: Pin<&mut ConnectedComponents>) -> Result<()>;
        fn hasFinished(self: &ConnectedComponents) -> bool;
        fn numberOfComponents(self: &ConnectedComponents) -> u64;
        fn ConnectedComponentsGetComponentSizes(
            algo: &ConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn ConnectedComponentsGetComponents(
            algo: &ConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn ConnectedComponentsExtractLargestConnectedComponent(
            g: &Graph,
            compact_graph: bool,
        ) -> UniquePtr<Graph>;
        fn componentOfNode(self: &ConnectedComponents, u: u64) -> u64;
        fn ConnectedComponentsGetPartition(algo: &ConnectedComponents) -> UniquePtr<Partition>;

        type DynConnectedComponents;
        fn NewDynConnectedComponents(g: &Graph) -> UniquePtr<DynConnectedComponents>;
        fn run(self: Pin<&mut DynConnectedComponents>) -> Result<()>;
        fn hasFinished(self: &DynConnectedComponents) -> bool;
        fn numberOfComponents(self: &DynConnectedComponents) -> u64;
        fn DynConnectedComponentsGetComponentSizes(
            algo: &DynConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn DynConnectedComponentsGetComponents(
            algo: &DynConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn componentOfNode(self: &DynConnectedComponents, u: u64) -> u64;
        fn DynConnectedComponentsGetPartition(
            algo: &DynConnectedComponents,
        ) -> UniquePtr<Partition>;
        fn DynConnectedComponentsUpdate(
            algo: Pin<&mut DynConnectedComponents>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        );
        fn DynConnectedComponentsUpdateBatch(
            algo: Pin<&mut DynConnectedComponents>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );

        type DynWeaklyConnectedComponents;
        fn NewDynWeaklyConnectedComponents(g: &Graph) -> UniquePtr<DynWeaklyConnectedComponents>;
        fn run(self: Pin<&mut DynWeaklyConnectedComponents>) -> Result<()>;
        fn hasFinished(self: &DynWeaklyConnectedComponents) -> bool;
        fn numberOfComponents(self: &DynWeaklyConnectedComponents) -> u64;
        fn DynWeaklyConnectedComponentsGetComponentSizes(
            algo: &DynWeaklyConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn DynWeaklyConnectedComponentsGetComponents(
            algo: &DynWeaklyConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn componentOfNode(self: &DynWeaklyConnectedComponents, u: u64) -> u64;
        fn DynWeaklyConnectedComponentsGetPartition(
            algo: &DynWeaklyConnectedComponents,
        ) -> UniquePtr<Partition>;
        fn DynWeaklyConnectedComponentsUpdate(
            algo: Pin<&mut DynWeaklyConnectedComponents>,
            kind: u8,
            u: u64,
            v: u64,
            ew: f64,
        );
        fn DynWeaklyConnectedComponentsUpdateBatch(
            algo: Pin<&mut DynWeaklyConnectedComponents>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );

        type ParallelConnectedComponents;
        fn NewParallelConnectedComponents(
            g: &Graph,
            coarsening: bool,
        ) -> UniquePtr<ParallelConnectedComponents>;
        fn run(self: Pin<&mut ParallelConnectedComponents>) -> Result<()>;
        fn hasFinished(self: &ParallelConnectedComponents) -> bool;
        fn numberOfComponents(self: &ParallelConnectedComponents) -> u64;
        fn ParallelConnectedComponentsGetComponentSizes(
            algo: &ParallelConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn ParallelConnectedComponentsGetComponents(
            algo: &ParallelConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn componentOfNode(self: &ParallelConnectedComponents, u: u64) -> u64;
        fn ParallelConnectedComponentsGetPartition(
            algo: &ParallelConnectedComponents,
        ) -> UniquePtr<Partition>;

        type StronglyConnectedComponents;
        fn NewStronglyConnectedComponents(g: &Graph) -> UniquePtr<StronglyConnectedComponents>;
        fn run(self: Pin<&mut StronglyConnectedComponents>) -> Result<()>;
        fn hasFinished(self: &StronglyConnectedComponents) -> bool;
        fn numberOfComponents(self: &StronglyConnectedComponents) -> u64;
        fn StronglyConnectedComponentsGetComponentSizes(
            algo: &StronglyConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn StronglyConnectedComponentsGetComponents(
            algo: &StronglyConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn componentOfNode(self: &StronglyConnectedComponents, u: u64) -> u64;
        fn StronglyConnectedComponentsGetPartition(
            algo: &StronglyConnectedComponents,
        ) -> UniquePtr<Partition>;

        type WeaklyConnectedComponents;
        fn NewWeaklyConnectedComponents(g: &Graph) -> UniquePtr<WeaklyConnectedComponents>;
        fn run(self: Pin<&mut WeaklyConnectedComponents>) -> Result<()>;
        fn hasFinished(self: &WeaklyConnectedComponents) -> bool;
        fn numberOfComponents(self: &WeaklyConnectedComponents) -> u64;
        fn WeaklyConnectedComponentsGetComponentSizes(
            algo: &WeaklyConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn WeaklyConnectedComponentsGetComponents(
            algo: &WeaklyConnectedComponents,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );
        fn componentOfNode(self: &WeaklyConnectedComponents, u: u64) -> u64;
        fn WeaklyConnectedComponentsGetPartition(
            algo: &WeaklyConnectedComponents,
        ) -> UniquePtr<Partition>;

        // ---- CORRELATION ----

        type Assortativity;
        fn NewAssortativity(g: &Graph, attributes: &[f64]) -> UniquePtr<Assortativity>;
        fn run(self: Pin<&mut Assortativity>) -> Result<()>;
        fn hasFinished(self: &Assortativity) -> bool;
        fn getCoefficient(self: &Assortativity) -> f64;

        // ---- DISTANCE ----

        type APSP;
        fn NewAPSP(g: &Graph) -> UniquePtr<APSP>;
        fn run(self: Pin<&mut APSP>) -> Result<()>;
        fn hasFinished(self: &APSP) -> bool;
        fn getDistance(self: &APSP, u: u64, v: u64) -> f64;
        fn APSPGetDistances(algo: &APSP, wt: &mut Vec<f64>) -> u64;

        type AStar;
        fn NewAStar(
            g: &Graph,
            heu: &CxxVector<f64>,
            src: u64,
            dst: u64,
            store_pred: bool,
        ) -> UniquePtr<AStar>;
        fn run(self: Pin<&mut AStar>) -> Result<()>;
        fn hasFinished(self: &AStar) -> bool;
        fn AStarGetPath(algo: &AStar) -> UniquePtr<CxxVector<u64>>;
        fn AStarGetPredecessors(algo: &AStar) -> UniquePtr<CxxVector<u64>>;
        fn getDistance(self: &AStar) -> f64;
        fn AStarGetDistances(algo: &AStar) -> UniquePtr<CxxVector<f64>>;
        fn setSource(self: Pin<&mut AStar>, src: u64);
        fn setTarget(self: Pin<&mut AStar>, dst: u64);
        fn AStarSetTargets(algo: Pin<&mut AStar>, dst: &[u64]);
        fn AStarGetTargetIndexMap(algo: &AStar, ks: &mut Vec<u64>, vs: &mut Vec<u64>);

        type AdamicAdarDistance;
        fn NewAdamicAdarDistance(g: &Graph) -> UniquePtr<AdamicAdarDistance>;
        fn preprocess(self: Pin<&mut AdamicAdarDistance>);
        fn distance(self: Pin<&mut AdamicAdarDistance>, u: u64, v: u64) -> f64;
        fn AdamicAdarDistanceGetEdgeScores(algo: &AdamicAdarDistance) -> UniquePtr<CxxVector<f64>>;

        type AlgebraicDistance;
        fn NewAlgebraicDistance(
            g: &Graph,
            n_systems: u64,
            n_iterations: u64,
            omega: f64,
            norm: u64,
            with_edge_scores: bool,
        ) -> UniquePtr<AlgebraicDistance>;
        fn preprocess(self: Pin<&mut AlgebraicDistance>);
        fn distance(self: Pin<&mut AlgebraicDistance>, u: u64, v: u64) -> f64;
        fn AlgebraicDistanceGetEdgeScores(algo: &AlgebraicDistance) -> UniquePtr<CxxVector<f64>>;

        type AllSimplePaths;
        fn NewAllSimplePaths(
            g: &Graph,
            src: u64,
            dst: u64,
            cut_off: u64,
        ) -> UniquePtr<AllSimplePaths>;
        fn run(self: Pin<&mut AllSimplePaths>) -> Result<()>;
        fn hasFinished(self: &AllSimplePaths) -> bool;
        fn numberOfSimplePaths(self: Pin<&mut AllSimplePaths>) -> u64;
        fn AllSimplePathsGetAllSimplePaths(algo: Pin<&mut AllSimplePaths>, vs: &mut Vec<u64>);

        type BFS;
        fn NewBFS(
            g: &Graph,
            src: u64,
            store_paths: bool,
            store_nodes_sorted_by_distance: bool,
            dst: u64,
        ) -> UniquePtr<BFS>;
        fn run(self: Pin<&mut BFS>) -> Result<()>;
        fn hasFinished(self: &BFS) -> bool;
        fn distance(self: &BFS, t: u64) -> f64;
        fn BFSGetDistances(algo: Pin<&mut BFS>) -> UniquePtr<CxxVector<f64>>;
        fn _numberOfPaths(self: &BFS, t: u64) -> f64;
        fn BFSGetPredecessors(algo: &BFS, t: u64) -> UniquePtr<CxxVector<u64>>;
        fn BFSGetPath(algo: &BFS, t: u64, forward: bool) -> UniquePtr<CxxVector<u64>>;
        fn BFSGetPaths(algo: &BFS, t: u64, forward: bool, vs: &mut Vec<u64>);
        fn BFSGetNodeSortedByDistance(algo: &BFS) -> UniquePtr<CxxVector<u64>>;
        fn getReachableNodes(self: &BFS) -> u64;
        fn setSource(self: Pin<&mut BFS>, src: u64);
        fn setTarget(self: Pin<&mut BFS>, dst: u64);
        fn getSumOfDistances(self: &BFS) -> f64;

        type BidirectionalBFS;
        fn NewBidirectionalBFS(
            g: &Graph,
            src: u64,
            dst: u64,
            store_pred: bool,
        ) -> UniquePtr<BidirectionalBFS>;
        fn run(self: Pin<&mut BidirectionalBFS>) -> Result<()>;
        fn hasFinished(self: &BidirectionalBFS) -> bool;
        fn BidirectionalBFSGetPath(algo: &BidirectionalBFS) -> UniquePtr<CxxVector<u64>>;
        fn BidirectionalBFSGetPredecessors(algo: &BidirectionalBFS) -> UniquePtr<CxxVector<u64>>;
        fn getDistance(self: &BidirectionalBFS) -> f64;
        fn BidirectionalBFSGetDistances(algo: &BidirectionalBFS) -> UniquePtr<CxxVector<f64>>;
        fn setSource(self: Pin<&mut BidirectionalBFS>, src: u64);
        fn setTarget(self: Pin<&mut BidirectionalBFS>, dst: u64);
        fn BidirectionalBFSSetTargets(algo: Pin<&mut BidirectionalBFS>, dst: &[u64]);
        fn BidirectionalBFSGetTargetIndexMap(
            algo: &BidirectionalBFS,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );

        type BidirectionalDijkstra;
        fn NewBidirectionalDijkstra(
            g: &Graph,
            src: u64,
            dst: u64,
            store_pred: bool,
        ) -> UniquePtr<BidirectionalDijkstra>;
        fn run(self: Pin<&mut BidirectionalDijkstra>) -> Result<()>;
        fn hasFinished(self: &BidirectionalDijkstra) -> bool;
        fn BidirectionalDijkstraGetPath(algo: &BidirectionalDijkstra) -> UniquePtr<CxxVector<u64>>;
        fn BidirectionalDijkstraGetPredecessors(
            algo: &BidirectionalDijkstra,
        ) -> UniquePtr<CxxVector<u64>>;
        fn getDistance(self: &BidirectionalDijkstra) -> f64;
        fn BidirectionalDijkstraGetDistances(
            algo: &BidirectionalDijkstra,
        ) -> UniquePtr<CxxVector<f64>>;
        fn setSource(self: Pin<&mut BidirectionalDijkstra>, src: u64);
        fn setTarget(self: Pin<&mut BidirectionalDijkstra>, dst: u64);
        fn BidirectionalDijkstraSetTargets(algo: Pin<&mut BidirectionalDijkstra>, dst: &[u64]);
        fn BidirectionalDijkstraGetTargetIndexMap(
            algo: &BidirectionalDijkstra,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );

        type CommuteTimeDistance;
        fn NewCommuteTimeDistance(g: &Graph, tol: f64) -> UniquePtr<CommuteTimeDistance>;
        fn run(self: Pin<&mut CommuteTimeDistance>) -> Result<()>;
        fn hasFinished(self: &CommuteTimeDistance) -> bool;
        fn runApproximation(self: Pin<&mut CommuteTimeDistance>);
        fn runParallelApproximation(self: Pin<&mut CommuteTimeDistance>);
        fn getSetupTime(self: &CommuteTimeDistance) -> u64;
        fn runSinglePair(self: Pin<&mut CommuteTimeDistance>, u: u64, v: u64) -> f64;
        fn distance(self: Pin<&mut CommuteTimeDistance>, u: u64, v: u64) -> f64;
        fn runSingleSource(self: Pin<&mut CommuteTimeDistance>, u: u64) -> f64;

        type Diameter;
        fn NewDiameter(g: &Graph, algo_t: u8, error: f64, n_samples: u64) -> UniquePtr<Diameter>;
        fn run(self: Pin<&mut Diameter>) -> Result<()>;
        fn hasFinished(self: &Diameter) -> bool;
        fn DiameterGetDiameter(algo: &Diameter, lower: &mut u64, upper: &mut u64);

        type Dijkstra;
        fn NewDijkstra(
            g: &Graph,
            src: u64,
            store_paths: bool,
            store_nodes_sorted_by_distance: bool,
            dst: u64,
        ) -> UniquePtr<Dijkstra>;
        fn run(self: Pin<&mut Dijkstra>) -> Result<()>;
        fn hasFinished(self: &Dijkstra) -> bool;
        fn distance(self: &Dijkstra, t: u64) -> f64;
        fn DijkstraGetDistances(algo: Pin<&mut Dijkstra>) -> UniquePtr<CxxVector<f64>>;
        fn _numberOfPaths(self: &Dijkstra, t: u64) -> f64;
        fn DijkstraGetPredecessors(algo: &Dijkstra, t: u64) -> UniquePtr<CxxVector<u64>>;
        fn DijkstraGetPath(algo: &Dijkstra, t: u64, forward: bool) -> UniquePtr<CxxVector<u64>>;
        fn DijkstraGetPaths(algo: &Dijkstra, t: u64, forward: bool, vs: &mut Vec<u64>);
        fn DijkstraGetNodeSortedByDistance(algo: &Dijkstra) -> UniquePtr<CxxVector<u64>>;
        fn getReachableNodes(self: &Dijkstra) -> u64;
        fn setSource(self: Pin<&mut Dijkstra>, src: u64);
        fn setTarget(self: Pin<&mut Dijkstra>, dst: u64);
        fn getSumOfDistances(self: &Dijkstra) -> f64;

        type DynAPSP;
        fn NewDynAPSP(g: Pin<&mut Graph>) -> UniquePtr<DynAPSP>;
        fn run(self: Pin<&mut DynAPSP>) -> Result<()>;
        fn hasFinished(self: &DynAPSP) -> bool;
        fn getDistance(self: &DynAPSP, u: u64, v: u64) -> f64;
        fn DynAPSPGetDistances(algo: &DynAPSP, wt: &mut Vec<f64>) -> u64;
        fn DynAPSPUpdate(algo: Pin<&mut DynAPSP>, kind: u8, u: u64, v: u64, ew: f64);
        fn DynAPSPUpdateBatch(
            algo: Pin<&mut DynAPSP>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );

        type DynBFS;
        fn NewDynBFS(g: &Graph, src: u64, store_predecessors: bool) -> UniquePtr<DynBFS>;
        fn run(self: Pin<&mut DynBFS>) -> Result<()>;
        fn hasFinished(self: &DynBFS) -> bool;
        fn distance(self: &DynBFS, t: u64) -> f64;
        fn DynBFSGetDistances(algo: Pin<&mut DynBFS>) -> UniquePtr<CxxVector<f64>>;
        fn _numberOfPaths(self: &DynBFS, t: u64) -> f64;
        fn DynBFSGetPredecessors(algo: &DynBFS, t: u64) -> UniquePtr<CxxVector<u64>>;
        fn DynBFSGetPath(algo: &DynBFS, t: u64, forward: bool) -> UniquePtr<CxxVector<u64>>;
        fn DynBFSGetPaths(algo: &DynBFS, t: u64, forward: bool, vs: &mut Vec<u64>);
        fn DynBFSGetNodeSortedByDistance(algo: &DynBFS) -> UniquePtr<CxxVector<u64>>;
        fn getReachableNodes(self: &DynBFS) -> u64;
        fn setSource(self: Pin<&mut DynBFS>, src: u64);
        fn setTarget(self: Pin<&mut DynBFS>, dst: u64);
        fn getSumOfDistances(self: &DynBFS) -> f64;
        fn DynBFSUpdate(algo: Pin<&mut DynBFS>, kind: u8, u: u64, v: u64, ew: f64);
        fn DynBFSUpdateBatch(
            algo: Pin<&mut DynBFS>,
            kinds: &[u8],
            us: &[u64],
            vs: &[u64],
            ew: &[f64],
        );
        fn modified(self: Pin<&mut DynBFS>) -> bool;
        fn setTargetNode(self: Pin<&mut DynBFS>, t: u64);

        fn EccentricityGetValue(g: &Graph, u: u64, fartherest: &mut u64, dist: &mut u64);

        type EffectiveDiameter;
        fn NewEffectiveDiameter(g: &Graph, ratio: f64) -> UniquePtr<EffectiveDiameter>;
        fn run(self: Pin<&mut EffectiveDiameter>) -> Result<()>;
        fn hasFinished(self: &EffectiveDiameter) -> bool;
        fn getEffectiveDiameter(self: &EffectiveDiameter) -> f64;

        type EffectiveDiameterApproximation;
        fn NewEffectiveDiameterApproximation(
            g: &Graph,
            ratio: f64,
            k: u64,
            r: u64,
        ) -> UniquePtr<EffectiveDiameterApproximation>;
        fn run(self: Pin<&mut EffectiveDiameterApproximation>) -> Result<()>;
        fn hasFinished(self: &EffectiveDiameterApproximation) -> bool;
        fn getEffectiveDiameter(self: &EffectiveDiameterApproximation) -> f64;

        type HopPlotApproximation;
        fn NewHopPlotApproximation(
            g: &Graph,
            max_distance: u64,
            k: u64,
            r: u64,
        ) -> UniquePtr<HopPlotApproximation>;
        fn run(self: Pin<&mut HopPlotApproximation>) -> Result<()>;
        fn hasFinished(self: &HopPlotApproximation) -> bool;
        fn HopPlotApproximationGetHopPlot(
            algo: &HopPlotApproximation,
            ks: &mut Vec<u64>,
            vs: &mut Vec<f64>,
        );

        type JaccardDistance;
        fn NewJaccardDistance(g: &Graph, triangles: &CxxVector<u64>) -> UniquePtr<JaccardDistance>;
        fn preprocess(self: Pin<&mut JaccardDistance>);
        fn distance(self: Pin<&mut JaccardDistance>, u: u64, v: u64) -> f64;
        fn JaccardDistanceGetEdgeScores(algo: &JaccardDistance) -> UniquePtr<CxxVector<f64>>;

        type MultiTargetBFS;
        fn NewMultiTargetBFS(g: &Graph, src: u64, targets: &[u64]) -> UniquePtr<MultiTargetBFS>;
        fn run(self: Pin<&mut MultiTargetBFS>) -> Result<()>;
        fn hasFinished(self: &MultiTargetBFS) -> bool;
        fn MultiTargetBFSGetPath(algo: &MultiTargetBFS) -> UniquePtr<CxxVector<u64>>;
        fn MultiTargetBFSGetPredecessors(algo: &MultiTargetBFS) -> UniquePtr<CxxVector<u64>>;
        fn getDistance(self: &MultiTargetBFS) -> f64;
        fn MultiTargetBFSGetDistances(algo: &MultiTargetBFS) -> UniquePtr<CxxVector<f64>>;
        fn setSource(self: Pin<&mut MultiTargetBFS>, src: u64);
        fn setTarget(self: Pin<&mut MultiTargetBFS>, dst: u64);
        fn MultiTargetBFSSetTargets(algo: Pin<&mut MultiTargetBFS>, dst: &[u64]);
        fn MultiTargetBFSGetTargetIndexMap(
            algo: &MultiTargetBFS,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );

        type MultiTargetDijkstra;
        fn NewMultiTargetDijkstra(
            g: &Graph,
            src: u64,
            targets: &[u64],
        ) -> UniquePtr<MultiTargetDijkstra>;
        fn run(self: Pin<&mut MultiTargetDijkstra>) -> Result<()>;
        fn hasFinished(self: &MultiTargetDijkstra) -> bool;
        fn MultiTargetDijkstraGetPath(algo: &MultiTargetDijkstra) -> UniquePtr<CxxVector<u64>>;
        fn MultiTargetDijkstraGetPredecessors(
            algo: &MultiTargetDijkstra,
        ) -> UniquePtr<CxxVector<u64>>;
        fn getDistance(self: &MultiTargetDijkstra) -> f64;
        fn MultiTargetDijkstraGetDistances(algo: &MultiTargetDijkstra)
            -> UniquePtr<CxxVector<f64>>;
        fn setSource(self: Pin<&mut MultiTargetDijkstra>, src: u64);
        fn setTarget(self: Pin<&mut MultiTargetDijkstra>, dst: u64);
        fn MultiTargetDijkstraSetTargets(algo: Pin<&mut MultiTargetDijkstra>, dst: &[u64]);
        fn MultiTargetDijkstraGetTargetIndexMap(
            algo: &MultiTargetDijkstra,
            ks: &mut Vec<u64>,
            vs: &mut Vec<u64>,
        );

        type NeighborhoodFunction;
        fn NewNeighborhoodFunction(g: &Graph) -> UniquePtr<NeighborhoodFunction>;
        fn run(self: Pin<&mut NeighborhoodFunction>) -> Result<()>;
        fn hasFinished(self: &NeighborhoodFunction) -> bool;
        fn NeighborhoodFunctionGetNeighborhoodFunction(
            algo: &NeighborhoodFunction,
        ) -> UniquePtr<CxxVector<u64>>;

        type NeighborhoodFunctionApproximation;
        fn NewNeighborhoodFunctionApproximation(
            g: &Graph,
            k: u64,
            r: u64,
        ) -> UniquePtr<NeighborhoodFunctionApproximation>;
        fn run(self: Pin<&mut NeighborhoodFunctionApproximation>) -> Result<()>;
        fn hasFinished(self: &NeighborhoodFunctionApproximation) -> bool;
        fn NeighborhoodFunctionApproximationGetNeighborhoodFunction(
            algo: &NeighborhoodFunctionApproximation,
        ) -> UniquePtr<CxxVector<u64>>;

        type NeighborhoodFunctionHeuristic;
        fn NewNeighborhoodFunctionHeuristic(
            g: &Graph,
            n_samples: u64,
            strategy: u8,
        ) -> UniquePtr<NeighborhoodFunctionHeuristic>;
        fn run(self: Pin<&mut NeighborhoodFunctionHeuristic>) -> Result<()>;
        fn hasFinished(self: &NeighborhoodFunctionHeuristic) -> bool;
        fn NeighborhoodFunctionHeuristicGetNeighborhoodFunction(
            algo: &NeighborhoodFunctionHeuristic,
        ) -> UniquePtr<CxxVector<u64>>;

        type PrunedLandmarkLabeling;
        fn NewPrunedLandmarkLabeling(g: &Graph) -> UniquePtr<PrunedLandmarkLabeling>;
        fn run(self: Pin<&mut PrunedLandmarkLabeling>) -> Result<()>;
        fn hasFinished(self: &PrunedLandmarkLabeling) -> bool;
        fn query(self: &PrunedLandmarkLabeling, u: u64, v: u64) -> u64;

        type ReverseBFS;
        fn NewReverseBFS(
            g: &Graph,
            src: u64,
            store_paths: bool,
            store_nodes_sorted_by_distance: bool,
            dst: u64,
        ) -> UniquePtr<ReverseBFS>;
        fn run(self: Pin<&mut ReverseBFS>) -> Result<()>;
        fn hasFinished(self: &ReverseBFS) -> bool;
        fn distance(self: &ReverseBFS, t: u64) -> f64;
        fn ReverseBFSGetDistances(algo: Pin<&mut ReverseBFS>) -> UniquePtr<CxxVector<f64>>;
        fn _numberOfPaths(self: &ReverseBFS, t: u64) -> f64;
        fn ReverseBFSGetPredecessors(algo: &ReverseBFS, t: u64) -> UniquePtr<CxxVector<u64>>;
        fn ReverseBFSGetPath(algo: &ReverseBFS, t: u64, forward: bool)
            -> UniquePtr<CxxVector<u64>>;
        fn ReverseBFSGetPaths(algo: &ReverseBFS, t: u64, forward: bool, vs: &mut Vec<u64>);
        fn ReverseBFSGetNodeSortedByDistance(algo: &ReverseBFS) -> UniquePtr<CxxVector<u64>>;
        fn getReachableNodes(self: &ReverseBFS) -> u64;
        fn setSource(self: Pin<&mut ReverseBFS>, src: u64);
        fn setTarget(self: Pin<&mut ReverseBFS>, dst: u64);
        fn getSumOfDistances(self: &ReverseBFS) -> f64;

        fn VolumeVolume(g: &Graph, r: f64, n_samples: u64) -> f64;
        fn VolumeVolumes(g: &Graph, rs: &[f64], n_samples: u64) -> UniquePtr<CxxVector<f64>>;

        // ---- EMBEDDING ----

        type Node2Vec;
        fn NewNode2Vec(g: &Graph, p: f64, q: f64, l: u64, n: u64, d: u64) -> UniquePtr<Node2Vec>;
        fn run(self: Pin<&mut Node2Vec>) -> Result<()>;
        fn hasFinished(self: &Node2Vec) -> bool;
        fn Node2VecGetFeatures(algo: &Node2Vec, ret: &mut Vec<f64>);

        // ---- FLOW ----

        type EdmondsKarp;
        fn NewEdmondsKarp(g: &Graph, source: u64, sink: u64) -> UniquePtr<EdmondsKarp>;
        fn run(self: Pin<&mut EdmondsKarp>) -> Result<()>;
        fn hasFinished(self: &EdmondsKarp) -> bool;
        fn getMaxFlow(self: &EdmondsKarp) -> f64;
        fn EdmondsKarpGetSourceSet(algo: &EdmondsKarp) -> UniquePtr<CxxVector<u64>>;
        fn getFlow(self: &EdmondsKarp, u: u64, v: u64) -> f64;
        #[rust_name = "getEdgeFlow"]
        fn getFlow(self: &EdmondsKarp, e: u64) -> f64;
        fn EdmondsKarpGetFlowVector(algo: &EdmondsKarp) -> UniquePtr<CxxVector<f64>>;

        // ---- GENERATORS ----

        type BarabasiAlbertGenerator;
        fn NewBarabasiAlbertGenerator(
            k: u64,
            n_max: u64,
            n0: u64,
            batagelj: bool,
        ) -> UniquePtr<BarabasiAlbertGenerator>;
        fn BarabasiAlbertGeneratorGenerate(
            algo: Pin<&mut BarabasiAlbertGenerator>,
        ) -> UniquePtr<Graph>;

        type ChungLuGenerator;
        fn NewChungLuGenerator(degree_sequence: &CxxVector<u64>) -> UniquePtr<ChungLuGenerator>;
        fn ChungLuGeneratorGenerate(algo: Pin<&mut ChungLuGenerator>) -> UniquePtr<Graph>;

        type ClusteredRandomGraphGenerator;
        fn NewClusteredRandomGraphGenerator(
            n: u64,
            k: u64,
            p_intra: f64,
            p_inter: f64,
        ) -> UniquePtr<ClusteredRandomGraphGenerator>;
        fn ClusteredRandomGraphGeneratorGenerate(
            algo: Pin<&mut ClusteredRandomGraphGenerator>,
        ) -> UniquePtr<Graph>;
        fn ClusteredRandomGraphGeneratorGetCommunities(
            algo: Pin<&mut ClusteredRandomGraphGenerator>,
        ) -> UniquePtr<Partition>;

        type DorogovtsevMendesGenerator;
        fn NewDorogovtsevMendesGenerator(n_nodes: u64) -> UniquePtr<DorogovtsevMendesGenerator>;
        fn DorogovtsevMendesGeneratorGenerate(
            algo: Pin<&mut DorogovtsevMendesGenerator>,
        ) -> UniquePtr<Graph>;

        type DynamicDorogovtsevMendesGenerator;
        fn NewDynamicDorogovtsevMendesGenerator() -> UniquePtr<DynamicDorogovtsevMendesGenerator>;
        fn DynamicDorogovtsevMendesGeneratorGenerate(
            algo: Pin<&mut DynamicDorogovtsevMendesGenerator>,
            n_steps: u64,
            tps: &mut Vec<u8>,
            us: &mut Vec<u64>,
            vs: &mut Vec<u64>,
            ws: &mut Vec<f64>,
        );

        type DynamicForestFireGenerator;
        fn NewDynamicForestFireGenerator(
            p: f64,
            directed: bool,
            r: f64,
        ) -> UniquePtr<DynamicForestFireGenerator>;
        fn DynamicForestFireGeneratorGenerate(
            algo: Pin<&mut DynamicForestFireGenerator>,
            n_steps: u64,
            tps: &mut Vec<u8>,
            us: &mut Vec<u64>,
            vs: &mut Vec<u64>,
            ws: &mut Vec<f64>,
        );

        type DynamicHyperbolicGenerator;
        fn NewDynamicHyperbolicGenerator(
            n: u64,
            avg_degree: f64,
            exp: f64,
            t: f64,
            move_each_step: f64,
            move_distance: f64,
        ) -> UniquePtr<DynamicHyperbolicGenerator>;
        fn DynamicHyperbolicGeneratorGenerate(
            algo: Pin<&mut DynamicHyperbolicGenerator>,
            n_steps: u64,
            tps: &mut Vec<u8>,
            us: &mut Vec<u64>,
            vs: &mut Vec<u64>,
            ws: &mut Vec<f64>,
        );
        fn DynamicHyperbolicGeneratorGetGraph(
            algo: &DynamicHyperbolicGenerator,
        ) -> UniquePtr<Graph>;
        fn DynamicHyperbolicGeneratorGetCoordinates(
            algo: &DynamicHyperbolicGenerator,
            xs: &mut Vec<f64>,
            ys: &mut Vec<f64>,
        );

        type DynamicPathGenerator;
        fn NewDynamicPathGenerator() -> UniquePtr<DynamicPathGenerator>;
        fn DynamicPathGeneratorGenerate(
            algo: Pin<&mut DynamicPathGenerator>,
            n_steps: u64,
            tps: &mut Vec<u8>,
            us: &mut Vec<u64>,
            vs: &mut Vec<u64>,
            ws: &mut Vec<f64>,
        );

        type DynamicPubWebGenerator;
        fn NewDynamicPubWebGenerator(
            num_nodes: u64,
            num_dense_areas: u64,
            neighbourhood_radius: f64,
            max_num_neighbours: u64,
            write_initial_graph_to_stream: bool,
        ) -> UniquePtr<DynamicPubWebGenerator>;
        fn DynamicPubWebGeneratorGenerate(
            algo: Pin<&mut DynamicPubWebGenerator>,
            n_steps: u64,
            tps: &mut Vec<u8>,
            us: &mut Vec<u64>,
            vs: &mut Vec<u64>,
            ws: &mut Vec<f64>,
        );
        fn DynamicPubWebGeneratorGetGraph(algo: &DynamicPubWebGenerator) -> UniquePtr<Graph>;
        fn DynamicPubWebGeneratorGetCoordinates(
            algo: &DynamicPubWebGenerator,
            xs: &mut Vec<f64>,
            ys: &mut Vec<f64>,
        );
        fn DynamicPubWebGeneratorGetNewCoordinates(
            algo: &DynamicPubWebGenerator,
            ns: &mut Vec<u64>,
            xs: &mut Vec<f64>,
            ys: &mut Vec<f64>,
        );

        type EdgeSwitchingMarkovChainGenerator;
        fn NewEdgeSwitchingMarkovChainGenerator(
            sequence: &CxxVector<u64>,
            ignore_if_not_realizable: bool,
            num_switches_per_edge: u64,
        ) -> UniquePtr<EdgeSwitchingMarkovChainGenerator>;
        fn EdgeSwitchingMarkovChainGeneratorGenerate(
            algo: Pin<&mut EdgeSwitchingMarkovChainGenerator>,
        ) -> UniquePtr<Graph>;
        fn isRealizable(self: Pin<&mut EdgeSwitchingMarkovChainGenerator>) -> bool;
        fn getRealizable(self: &EdgeSwitchingMarkovChainGenerator) -> bool;

        type ErdosRenyiGenerator;
        fn NewErdosRenyiGenerator(
            n_nodes: u64,
            prob: f64,
            directed: bool,
            self_loops: bool,
        ) -> UniquePtr<ErdosRenyiGenerator>;
        fn ErdosRenyiGeneratorGenerate(algo: Pin<&mut ErdosRenyiGenerator>) -> UniquePtr<Graph>;

        type HavelHakimiGenerator;
        fn NewHavelHakimiGenerator(
            sequence: &CxxVector<u64>,
            ignore_if_not_realizable: bool,
        ) -> UniquePtr<HavelHakimiGenerator>;
        fn HavelHakimiGeneratorGenerate(algo: Pin<&mut HavelHakimiGenerator>) -> UniquePtr<Graph>;
        fn isRealizable(self: Pin<&mut HavelHakimiGenerator>) -> bool;
        fn getRealizable(self: &HavelHakimiGenerator) -> bool;

        type HyperbolicGenerator;
        fn NewHyperbolicGenerator(
            n: u64,
            avg_degree: f64,
            exp: f64,
            T: f64,
        ) -> UniquePtr<HyperbolicGenerator>;
        fn HyperbolicGeneratorGenerate(algo: Pin<&mut HyperbolicGenerator>) -> UniquePtr<Graph>;
        #[rust_name = "HyperbolicGeneratorGenerateAdvanced"]
        fn HyperbolicGeneratorGenerate(
            algo: Pin<&mut HyperbolicGenerator>,
            angles: &[f64],
            radii: &[f64],
            r: f64,
            T: f64,
        ) -> UniquePtr<Graph>;
        fn setLeafCapacity(self: Pin<&mut HyperbolicGenerator>, capacity: u64);
        fn setTheoreticalSplit(self: Pin<&mut HyperbolicGenerator>, split: bool);
        fn setBalance(self: Pin<&mut HyperbolicGenerator>, balance: f64);

        type LFRGenerator;
        fn NewLFRGenerator(n: u64) -> UniquePtr<LFRGenerator>;
        fn LFRGeneratorGenerate(algo: Pin<&mut LFRGenerator>) -> UniquePtr<Graph>;
        fn run(self: Pin<&mut LFRGenerator>) -> Result<()>;
        fn hasFinished(self: &LFRGenerator) -> bool;
        fn LFRGeneratorSetDegreeSequence(algo: Pin<&mut LFRGenerator>, seq: &[u64]);
        fn generatePowerlawDegreeSequence(
            self: Pin<&mut LFRGenerator>,
            avg_degree: u64,
            max_degree: u64,
            node_degree_exp: f64,
        );
        fn LFRGeneratorSetCommunitySizeSequence(algo: Pin<&mut LFRGenerator>, seq: &[u64]);
        fn LFRGeneratorSetPartition(algo: Pin<&mut LFRGenerator>, p: UniquePtr<Partition>);
        fn generatePowerlawCommunitySizeSequence(
            self: Pin<&mut LFRGenerator>,
            min_community_size: u64,
            max_community_size: u64,
            community_size_exp: f64,
        );
        fn setMu(self: Pin<&mut LFRGenerator>, mu: f64);
        #[rust_name = "setMuArray"]
        fn setMu(self: Pin<&mut LFRGenerator>, mu: &CxxVector<f64>);
        fn setMuWithBinomialDistribution(self: Pin<&mut LFRGenerator>, mu: f64);
        fn LFRGeneratorGetGraph(algo: &LFRGenerator) -> UniquePtr<Graph>;
        fn LFRGeneratorGetPartition(algo: &LFRGenerator) -> UniquePtr<Partition>;

        type MocnikGenerator;
        fn NewMocnikGenerator(
            dim: u64,
            n: u64,
            k: f64,
            weighted: bool,
        ) -> UniquePtr<MocnikGenerator>;
        fn MocnikGeneratorGenerate(algo: Pin<&mut MocnikGenerator>) -> UniquePtr<Graph>;

        type MocnikGeneratorBasic;
        fn NewMocnikGeneratorBasic(dim: u64, n: u64, k: f64) -> UniquePtr<MocnikGeneratorBasic>;
        fn MocnikGeneratorBasicGenerate(algo: Pin<&mut MocnikGeneratorBasic>) -> UniquePtr<Graph>;

        type PowerlawDegreeSequence;
        fn NewPowerlawDegreeSequence(
            min_deg: u64,
            max_deg: u64,
            gamma: f64,
        ) -> UniquePtr<PowerlawDegreeSequence>;
        fn run(self: Pin<&mut PowerlawDegreeSequence>) -> Result<()>;
        fn hasFinished(self: &PowerlawDegreeSequence) -> bool;
        fn PowerlawDegreeSequenceGetDegreeSequence(
            algo: &PowerlawDegreeSequence,
            num_nodes: u64,
        ) -> UniquePtr<CxxVector<u64>>;
        fn setMinimumFromAverageDegree(self: Pin<&mut PowerlawDegreeSequence>, avg_deg: f64);
        fn setGammaFromAverageDegree(
            self: Pin<&mut PowerlawDegreeSequence>,
            avg_deg: f64,
            min_gamma: f64,
            max_gamma: f64,
        );
        fn setMinimumDegree(self: Pin<&mut PowerlawDegreeSequence>, min_deg: u64);
        fn getMinimumDegree(self: &PowerlawDegreeSequence) -> u64;
        fn getMaximumDegree(self: &PowerlawDegreeSequence) -> u64;
        fn setGamma(self: Pin<&mut PowerlawDegreeSequence>, gamma: f64);
        fn getGamma(self: &PowerlawDegreeSequence) -> f64;
        fn getExpectedAverageDegree(self: &PowerlawDegreeSequence) -> f64;
        fn getDegree(self: &PowerlawDegreeSequence) -> u64;

        type PubWebGenerator;
        fn NewPubWebGenerator(
            num_nodes: u64,
            num_dense_areas: u64,
            neighbourhood_radius: f64,
            max_num_neighbours: u64,
        ) -> UniquePtr<PubWebGenerator>;
        fn PubWebGeneratorGenerate(algo: Pin<&mut PubWebGenerator>) -> UniquePtr<Graph>;
        fn PubWebGeneratorGetCoordinates(
            algo: &PubWebGenerator,
            xs: &mut Vec<f64>,
            ys: &mut Vec<f64>,
        );

        type RegularRingLatticeGenerator;
        fn NewRegularRingLatticeGenerator(
            num_nodes: u64,
            num_neighbours: u64,
        ) -> UniquePtr<RegularRingLatticeGenerator>;
        fn RegularRingLatticeGeneratorGenerate(
            algo: Pin<&mut RegularRingLatticeGenerator>,
        ) -> UniquePtr<Graph>;

        type RmatGenerator;
        fn NewRmatGenerator(
            scale: u64,
            edge_factor: u64,
            a: f64,
            b: f64,
            c: f64,
            d: f64,
            weighted: bool,
            reduce_nodes: u64,
        ) -> UniquePtr<RmatGenerator>;
        fn RmatGeneratorGenerate(algo: Pin<&mut RmatGenerator>) -> UniquePtr<Graph>;

        type WattsStrogatzGenerator;
        fn NewWattsStrogatzGenerator(
            n_nodes: u64,
            n_neighbours: u64,
            p: f64,
        ) -> UniquePtr<WattsStrogatzGenerator>;
        fn WattsStrogatzGeneratorGenerate(
            algo: Pin<&mut WattsStrogatzGenerator>,
        ) -> UniquePtr<Graph>;

        // ---- INDEPENDENT SET ----
        type Luby;
        fn NewLuby() -> UniquePtr<Luby>;
        fn LubyRun(algo: Pin<&mut Luby>, g: &Graph, ret: &mut Vec<bool>);
        fn LubyIsIndependentSet(algo: &Luby, set: &[bool], g: &Graph) -> bool;

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
