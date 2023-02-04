use std::collections::BTreeMap;

use cxx::UniquePtr;
use miette::IntoDiagnostic;

use crate::{
    base::Algorithm,
    bridge::{self, *},
};

pub struct ApproximatePageRank {
    inner: UniquePtr<bridge::ApproximatePageRank>,
}

impl ApproximatePageRank {
    pub fn new(g: &crate::Graph, alpha: f64, epsilon: Option<f64>) -> Self {
        Self {
            inner: NewApproximatePageRank(g, alpha, epsilon.unwrap_or(1e-12)),
        }
    }
    pub fn run(&mut self, seeds: &[u64]) -> impl Iterator<Item = (u64, f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        ApproximatePageRankRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        ks.into_iter().zip(vs.into_iter())
    }
}

pub struct SelectiveCommunityDetectorBase {
    pub(crate) inner: UniquePtr<bridge::SelectiveCommunityDetector>,
}

pub trait SelectiveCommunityDetector {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>>;
    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64>;
    fn into_base(self) -> SelectiveCommunityDetectorBase;
}

pub struct CliqueDetect {
    inner: UniquePtr<bridge::CliqueDetect>,
}

impl CliqueDetect {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewCliqueDetect(g),
        }
    }
}

impl SelectiveCommunityDetector for CliqueDetect {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        CliqueDetectRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        CliqueDetectExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: CliqueDetectAsBase(self.inner),
        }
    }
}

pub struct CombinedSCD {
    inner: UniquePtr<bridge::CombinedSCD>,
}

impl CombinedSCD {
    pub fn new(
        g: &crate::Graph,
        left: &mut SelectiveCommunityDetectorBase,
        right: &mut SelectiveCommunityDetectorBase,
    ) -> Self {
        Self {
            inner: NewCombinedSCD(g, left.inner.pin_mut(), right.inner.pin_mut()),
        }
    }
}

impl SelectiveCommunityDetector for CombinedSCD {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        CombinedSCDRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        CombinedSCDExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: CombinedSCDAsBase(self.inner),
        }
    }
}

pub struct GCE {
    inner: UniquePtr<bridge::GCE>,
}

impl GCE {
    pub const M: &str = "M";
    pub const L: &str = "L";

    pub fn new(g: &crate::Graph, q: &str) -> Self {
        Self {
            inner: NewGCE(g, q),
        }
    }
}

impl SelectiveCommunityDetector for GCE {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        GCERun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        GCEExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: GCEAsBase(self.inner),
        }
    }
}

pub struct LFMLocal {
    inner: UniquePtr<bridge::LFMLocal>,
}

impl LFMLocal {
    pub fn new(g: &crate::Graph, alpha: Option<f64>) -> Self {
        Self {
            inner: NewLFMLocal(g, alpha.unwrap_or(1.)),
        }
    }
}

impl SelectiveCommunityDetector for LFMLocal {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        LFMLocalRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        LFMLocalExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: LFMLocalAsBase(self.inner),
        }
    }
}

pub struct LocalT {
    inner: UniquePtr<bridge::LocalT>,
}

impl LocalT {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewLocalT(g),
        }
    }
}

impl SelectiveCommunityDetector for LocalT {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        LocalTRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        LocalTExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: LocalTAsBase(self.inner),
        }
    }
}

pub struct LocalTightnessExpansion {
    inner: UniquePtr<bridge::LocalTightnessExpansion>,
}

impl LocalTightnessExpansion {
    pub fn new(g: &crate::Graph, alpha: Option<f64>) -> Self {
        Self {
            inner: NewLocalTightnessExpansion(g, alpha.unwrap_or(1.)),
        }
    }
}

impl SelectiveCommunityDetector for LocalTightnessExpansion {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        LocalTightnessExpansionRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        LocalTightnessExpansionExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: LocalTightnessExpansionAsBase(self.inner),
        }
    }
}

pub struct PageRankNibble {
    inner: UniquePtr<bridge::PageRankNibble>,
}

impl PageRankNibble {
    pub fn new(g: &crate::Graph, alpha: f64, epsilon: f64) -> Self {
        Self {
            inner: NewPageRankNibble(g, alpha, epsilon),
        }
    }
}

impl SelectiveCommunityDetector for PageRankNibble {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        PageRankNibbleRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        PageRankNibbleExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: PageRankNibbleAsBase(self.inner),
        }
    }
}

pub struct RandomBFS {
    inner: UniquePtr<bridge::RandomBFS>,
}

impl RandomBFS {
    pub fn new(g: &crate::Graph, c: &crate::Cover) -> Self {
        Self {
            inner: NewRandomBFS(g, c),
        }
    }
}

impl SelectiveCommunityDetector for RandomBFS {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        RandomBFSRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        RandomBFSExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: RandomBFSAsBase(self.inner),
        }
    }
}

pub struct SCDGroundTruthComparison {
    inner: UniquePtr<bridge::SCDGroundTruthComparison>,
}

impl SCDGroundTruthComparison {
    pub fn new(
        g: &crate::Graph,
        ground_truth: &crate::Cover,
        found: &BTreeMap<u64, Vec<u64>>,
        ignore_seeds: bool,
    ) -> Self {
        let mut ks = vec![];
        let mut vs = vec![];

        for (k, v) in found {
            for vv in v {
                ks.push(*k);
                vs.push(*vv);
            }
        }

        Self {
            inner: NewSCDGroundTruthComparison(g, ground_truth, &ks, &vs, ignore_seeds),
        }
    }

    pub fn get_individual_jaccard(&self) -> impl Iterator<Item = (u64, f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        SCDGroundTruthComparisonGetIndividualJaccard(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs.into_iter())
    }

    pub fn get_individual_precision(&self) -> impl Iterator<Item = (u64, f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        SCDGroundTruthComparisonGetIndividualPrecision(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs.into_iter())
    }

    pub fn get_individual_recall(&self) -> impl Iterator<Item = (u64, f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        SCDGroundTruthComparisonGetIndividualRecall(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs.into_iter())
    }

    pub fn get_individual_f1(&self) -> impl Iterator<Item = (u64, f64)> {
        let mut ks = vec![];
        let mut vs = vec![];
        SCDGroundTruthComparisonGetIndividualF1(&self.inner, &mut ks, &mut vs);
        ks.into_iter().zip(vs.into_iter())
    }

    pub fn get_average_jaccard(&self) -> f64 {
        self.inner.getAverageJaccard()
    }

    pub fn get_average_precision(&self) -> f64 {
        self.inner.getAveragePrecision()
    }

    pub fn get_average_recall(&self) -> f64 {
        self.inner.getAverageRecall()
    }

    pub fn get_average_f1(&self) -> f64 {
        self.inner.getAverageF1()
    }
}

impl Algorithm for SCDGroundTruthComparison {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct SetConductance {
    inner: UniquePtr<bridge::SetConductance>,
}

impl SetConductance {
    pub fn new(g: &crate::Graph, community: &[u64]) -> Self {
        Self {
            inner: NewSetConductance(g, community),
        }
    }
    pub fn get_conductance(&self) -> f64 {
        self.inner.getConductance()
    }
}

impl Algorithm for SetConductance {
    fn run(&mut self) -> miette::Result<()> {
        self.inner.pin_mut().run().into_diagnostic()
    }

    fn has_finished(&self) -> bool {
        self.inner.hasFinished()
    }
}

pub struct TCE {
    inner: UniquePtr<bridge::TCE>,
}

impl TCE {
    pub fn new(g: &crate::Graph, refine: bool, use_jaccard: bool) -> Self {
        Self {
            inner: NewTCE(g, refine, use_jaccard),
        }
    }
}

impl SelectiveCommunityDetector for TCE {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        TCERun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        TCEExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: TCEAsBase(self.inner),
        }
    }
}

pub struct TwoPhaseL {
    inner: UniquePtr<bridge::TwoPhaseL>,
}

impl TwoPhaseL {
    pub fn new(g: &crate::Graph) -> Self {
        Self {
            inner: NewTwoPhaseL(g),
        }
    }
}

impl SelectiveCommunityDetector for TwoPhaseL {
    fn run(&mut self, seeds: &[u64]) -> BTreeMap<u64, Vec<u64>> {
        let mut ks = vec![];
        let mut vs = vec![];
        TwoPhaseLRun(self.inner.pin_mut(), seeds, &mut ks, &mut vs);
        let mut ret: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (k, v) in ks.into_iter().zip(vs.into_iter()) {
            ret.entry(k).or_default().push(v);
        }
        ret
    }

    fn expand_one_community(&mut self, seeds: &[u64]) -> Vec<u64> {
        let mut ret = vec![];
        TwoPhaseLExpandOneCommunity(self.inner.pin_mut(), seeds, &mut ret);
        ret
    }

    fn into_base(self) -> SelectiveCommunityDetectorBase {
        SelectiveCommunityDetectorBase {
            inner: TwoPhaseLAsBase(self.inner),
        }
    }
}
