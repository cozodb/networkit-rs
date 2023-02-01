use std::ops::Deref;

use cxx::UniquePtr;

use crate::{
    bridge::{self, *},
    tools::NodeIter,
};

pub struct Cover {
    inner: UniquePtr<bridge::Cover>,
}

unsafe impl Send for Cover {}

impl Deref for Cover {
    type Target = bridge::Cover;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl Clone for Cover {
    fn clone(&self) -> Self {
        Self {
            inner: CopyCover(&self.inner),
        }
    }
}

impl From<&crate::Partition> for Cover {
    fn from(value: &crate::Partition) -> Self {
        Self {
            inner: NewCoverFromPartition(&value.inner),
        }
    }
}

impl Default for Cover {
    fn default() -> Self {
        Self { inner: NewCover() }
    }
}

impl Cover {
    pub fn new(z: u64) -> Self {
        Self {
            inner: NewCoverWithSize(z),
        }
    }
    pub fn add_to_subset(&mut self, s: u64, e: u64) {
        self.inner.pin_mut().addToSubset(s, e)
    }
    pub fn all_to_singletons(&mut self) {
        self.inner.pin_mut().allToSingletons()
    }
    pub fn contains(&self, e: u64) -> bool {
        self.inner.contains(e)
    }
    pub fn extend(&mut self) -> u64 {
        self.inner.pin_mut().extend()
    }
    pub fn get_members(&self, s: u64) -> Vec<u64> {
        let mut ret = vec![];
        CVGetMembers(self, s, &mut ret);
        ret
    }
    pub fn get_subset_ids(&self) -> Vec<u64> {
        let mut ret = vec![];
        CVGetSubsetIds(self, &mut ret);
        ret
    }
    pub fn in_same_subset(&self, e1: u64, e2: u64) -> bool {
        self.inner.inSameSubset(e1, e2)
    }
    pub fn lower_bound(&self) -> u64 {
        self.inner.lowerBound()
    }
    pub fn merge_subsets(&mut self, s: u64, t: u64) {
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
    pub fn remove_from_subset(&mut self, s: u64, e: u64) {
        self.inner.pin_mut().removeFromSubset(s, e)
    }
    pub fn set_upper_bound(&mut self, upper: u64) {
        self.inner.pin_mut().setUpperBound(upper)
    }
    pub fn subset_size_map(&self) -> impl Iterator<Item = (u64, u64)> {
        let mut ks = vec![];
        let mut sz = vec![];
        CVSubsetSizeMap(self, &mut ks, &mut sz);
        ks.into_iter().zip(sz.into_iter())
    }
    pub fn subset_sizes(&self) -> impl Iterator<Item = u64> {
        NodeIter {
            nodes: CVSubsetSizes(self),
            at: 0,
        }
    }
    pub fn subsets_of(&self, e: u64) -> impl Iterator<Item = u64> {
        let ret = CVSubsetsOf(self, e);
        NodeIter { nodes: ret, at: 0 }
    }
    pub fn to_singleton(&mut self, e: u64) -> u64 {
        self.inner.pin_mut().toSingleton(e)
    }
    pub fn upper_bound(&self) -> u64 {
        self.inner.upperBound()
    }
}
