use std::ops::Deref;

use cxx::UniquePtr;

use crate::{
    bridge::{self, *},
    tools::NodeIter,
};

pub struct Partition {
    pub(crate) inner: UniquePtr<bridge::Partition>,
}

unsafe impl Send for Partition {}

impl Clone for Partition {
    fn clone(&self) -> Self {
        Self {
            inner: CopyPartition(self),
        }
    }
}

impl Deref for Partition {
    type Target = bridge::Partition;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl Default for Partition {
    fn default() -> Self {
        Partition::new(0)
    }
}

impl Partition {
    pub fn new(z: u64) -> Self {
        Self {
            inner: NewPartition(z),
        }
    }
    pub fn add_to_subset(&mut self, s: u64, e: u64) {
        self.inner.pin_mut().addToSubset(s, e)
    }
    pub fn all_to_singletons(&mut self) {
        self.inner.pin_mut().allToSingletons()
    }
    pub fn compact(&mut self, use_turbo: bool) {
        self.inner.pin_mut().compact(use_turbo)
    }
    pub fn contains(&self, e: u64) -> bool {
        self.inner.contains(e)
    }
    pub fn extend(&mut self) -> u64 {
        self.inner.pin_mut().extend()
    }
    pub fn get_members(&self, s: u64) -> Vec<u64> {
        let mut ret = vec![];
        PTGetMembers(self, s, &mut ret);
        ret
    }
    pub fn get_name(&self) -> String {
        let cs = PTGetName(self);
        cs.to_string()
    }
    pub fn get_subset_ids(&self) -> Vec<u64> {
        let mut ret = vec![];
        PTGetSubsetIds(self, &mut ret);
        ret
    }
    pub fn get_vector(&self) -> Vec<u64> {
        let r = self.inner.getVector();
        r.iter().cloned().collect()
    }
    pub fn in_same_subset(&self, e1: u64, e2: u64) -> bool {
        self.inner.inSameSubset(e1, e2)
    }
    pub fn lower_bound(&self) -> u64 {
        self.inner.lowerBound()
    }
    pub fn merge_subsets(&mut self, s: u64, t: u64) -> u64 {
        self.inner.pin_mut().mergeSubsets(s, t)
    }
    pub fn move_to_subset(&mut self, s: u64, e: u64) {
        self.inner.pin_mut().moveToSubset(s, e)
    }
    pub fn number_of_elements(&self) -> u64 {
        self.inner.numberOfElements()
    }
    pub fn number_of_subsets(&self) -> u64 {
        self.inner.numberOfSubsets()
    }
    pub fn set_name(&mut self, name: &str) {
        PTSetName(self.inner.pin_mut(), name)
    }
    pub fn set_upper_bound(&mut self, upper: u64) {
        self.inner.pin_mut().setUpperBound(upper)
    }
    pub fn subset_of(&self, e: u64) -> u64 {
        self.inner.subsetOf(e)
    }
    pub fn subset_size_map(&self) -> impl Iterator<Item = (u64, u64)> {
        let mut ks = vec![];
        let mut sz = vec![];
        PTSubsetSizeMap(self, &mut ks, &mut sz);
        ks.into_iter().zip(sz.into_iter())
    }
    pub fn subset_sizes(&self) -> impl Iterator<Item = u64> {
        NodeIter {
            nodes: PTSubsetSizes(self),
            at: 0,
        }
    }
    pub fn to_singleton(&mut self, e: u64) {
        self.inner.pin_mut().toSingleton(e)
    }
    pub fn upper_bound(&self) -> u64 {
        self.inner.upperBound()
    }
}
