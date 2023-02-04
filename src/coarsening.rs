use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
    tools::NodeIter,
};

pub trait GraphCoarsening: Algorithm {
    fn get_coarse_graph(&self) -> crate::Graph;
    fn get_fine_to_coarse_node_mapping(&self) -> NodeIter;
    fn get_coarse_to_fine_node_mapping(&self) -> BTreeMap<u64, Vec<u64>>;
}

pub struct ParallelPartitionCoarsening {
    inner: UniquePtr<bridge::ParallelPartitionCoarsening>,
}

impl ParallelPartitionCoarsening {
    pub fn new(g: &crate::Graph, zeta: &crate::Partition, parallel: bool) -> Self {
        Self {
            inner: NewParallelPartitionCoarsening(g, zeta, parallel),
        }
    }
}

impl GraphCoarsening for ParallelPartitionCoarsening {
    fn get_coarse_graph(&self) -> crate::Graph {
        ParallelPartitionCoarseningGetCoarseGraph(&self.inner).into()
    }

    fn get_fine_to_coarse_node_mapping(&self) -> NodeIter {
        NodeIter {
            nodes: ParallelPartitionCoarseningGetFineToCoarseNodeMapping(&self.inner),
            at: 0,
        }
    }

    fn get_coarse_to_fine_node_mapping(&self) -> BTreeMap<u64, Vec<u64>> {
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (i, k) in self.get_fine_to_coarse_node_mapping().enumerate() {
            ret.entry(k).or_default().push(i as u64);
        }
        ret
    }
}

impl Algorithm for ParallelPartitionCoarsening {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}
