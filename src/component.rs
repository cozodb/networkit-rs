use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::{Algorithm, DynAlgorithm},
    bridge::{self, *},
};

pub trait ComponentDecomposition {
    fn number_of_components(&self) -> u64;
    fn component_of_node(&self, u: u64) -> u64;
    fn get_partition(&self) -> crate::Partition;
    fn get_component_sizes(&self) -> BTreeMap<u64, u64>;
    fn get_components(&self) -> Vec<Vec<u64>>;
}



pub struct BiconnectedComponents {
    inner: UniquePtr<bridge::BiconnectedComponents>,
}

impl BiconnectedComponents {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewBiconnectedComponents(g),
        }
    }
    pub fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
    pub fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        BiconnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
    pub fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        BiconnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
    pub fn get_component_of_node(&self, u: u64) -> Vec<u64> {
        let mut vs = vec![];
        BiconnectedComponentsGetComponentOfNode(&self.inner,u, &mut vs);
        vs
    }
}

impl Algorithm for BiconnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}




pub struct ConnectedComponents {
    inner: UniquePtr<bridge::ConnectedComponents>,
}

impl ConnectedComponents {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewConnectedComponents(g),
        }
    }

    pub fn extract_largest_connected_component(g: &crate::Graph, compact_graph: bool) -> crate::Graph {
        ConnectedComponentsExtractLargestConnectedComponent(g, compact_graph).into()
    }
}

impl Algorithm for ConnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}


impl ComponentDecomposition  for ConnectedComponents{
    
     fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
     fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        ConnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
     fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        ConnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
     fn component_of_node(&self, u: u64) -> u64 {
        self.inner.componentOfNode(u)
    }

    fn get_partition(&self) -> crate::Partition {
        ConnectedComponentsGetPartition(&self.inner).into()
    }
}




pub struct DynConnectedComponents {
    inner: UniquePtr<bridge::DynConnectedComponents>,
}

impl DynConnectedComponents {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewDynConnectedComponents(g),
        }
    }
}

impl Algorithm for DynConnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}


impl ComponentDecomposition  for DynConnectedComponents{
    
     fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
     fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        DynConnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
     fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        DynConnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
     fn component_of_node(&self, u: u64) -> u64 {
        self.inner.componentOfNode(u)
    }

    fn get_partition(&self) -> crate::Partition {
        DynConnectedComponentsGetPartition(&self.inner).into()
    }
}



impl DynAlgorithm for DynConnectedComponents {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynConnectedComponentsUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynConnectedComponentsUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}




pub struct DynWeaklyConnectedComponents {
    inner: UniquePtr<bridge::DynWeaklyConnectedComponents>,
}

impl DynWeaklyConnectedComponents {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewDynWeaklyConnectedComponents(g),
        }
    }
}

impl Algorithm for DynWeaklyConnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}


impl ComponentDecomposition  for DynWeaklyConnectedComponents{
    
     fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
     fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        DynWeaklyConnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
     fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        DynWeaklyConnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
     fn component_of_node(&self, u: u64) -> u64 {
        self.inner.componentOfNode(u)
    }

    fn get_partition(&self) -> crate::Partition {
        DynWeaklyConnectedComponentsGetPartition(&self.inner).into()
    }
}



impl DynAlgorithm for DynWeaklyConnectedComponents {
    fn update(&mut self, e: crate::base::GraphEvent) {
        DynWeaklyConnectedComponentsUpdate(self.inner.pin_mut(), e.kind as u8, e.u, e.v, e.ew)
    }

    fn update_batch(&mut self, es: &[crate::base::GraphEvent]) {
        let mut kinds = Vec::with_capacity(es.len());
        let mut us = Vec::with_capacity(es.len());
        let mut vs = Vec::with_capacity(es.len());
        let mut ews = Vec::with_capacity(es.len());
        for ev in es {
            kinds.push(ev.kind as u8);
            us.push(ev.u);
            vs.push(ev.v);
            ews.push(ev.ew);
        }
        DynWeaklyConnectedComponentsUpdateBatch(self.inner.pin_mut(), &kinds, &us, &vs, &ews);
    }
}



pub struct ParallelConnectedComponents {
    inner: UniquePtr<bridge::ParallelConnectedComponents>,
}

impl ParallelConnectedComponents {
    pub fn new(g: &crate::Graph, coarsening: bool) -> Self {
        Self {
            inner: NewParallelConnectedComponents(g, coarsening),
        }
    }
}

impl Algorithm for ParallelConnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}


impl ComponentDecomposition  for ParallelConnectedComponents{
    
     fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
     fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        ParallelConnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
     fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        ParallelConnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
     fn component_of_node(&self, u: u64) -> u64 {
        self.inner.componentOfNode(u)
    }

    fn get_partition(&self) -> crate::Partition {
        ParallelConnectedComponentsGetPartition(&self.inner).into()
    }
}




pub struct StronglyConnectedComponents {
    inner: UniquePtr<bridge::StronglyConnectedComponents>,
}

impl StronglyConnectedComponents {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewStronglyConnectedComponents(g),
        }
    }
}

impl Algorithm for StronglyConnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}


impl ComponentDecomposition  for StronglyConnectedComponents{
    
     fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
     fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        StronglyConnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
     fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        StronglyConnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
     fn component_of_node(&self, u: u64) -> u64 {
        self.inner.componentOfNode(u)
    }

    fn get_partition(&self) -> crate::Partition {
        StronglyConnectedComponentsGetPartition(&self.inner).into()
    }
}




pub struct WeaklyConnectedComponents {
    inner: UniquePtr<bridge::WeaklyConnectedComponents>,
}

impl WeaklyConnectedComponents {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewWeaklyConnectedComponents(g),
        }
    }
}

impl Algorithm for WeaklyConnectedComponents {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}


impl ComponentDecomposition  for WeaklyConnectedComponents{
    
     fn number_of_components(&self) -> u64 {
        self.inner.numberOfComponents()
    }
     fn get_component_sizes(&self) -> BTreeMap<u64, u64> {
        let mut ks = vec![];
        let mut vs = vec![];
        WeaklyConnectedComponentsGetComponentSizes(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs).collect()
    }
     fn get_components(&self) -> Vec<Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        WeaklyConnectedComponentsGetComponents(&self.inner, &mut ks, &mut vs);
        if ks.is_empty() {
            vec![]
        } else {
            let l = ks.last().unwrap();
            let mut ret = vec![vec![]; *l as usize + 1];
            for (idx, v) in ks.into_iter().zip(vs) {
                ret[idx as usize].push(v);
            }
            ret
        }
    }
     fn component_of_node(&self, u: u64) -> u64 {
        self.inner.componentOfNode(u)
    }

    fn get_partition(&self) -> crate::Partition {
        WeaklyConnectedComponentsGetPartition(&self.inner).into()
    }
}
