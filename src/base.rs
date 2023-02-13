use miette::Result;

pub trait Algorithm {
    fn run(&mut self) -> Result<()>;
    fn has_finished(&self) -> bool;
}

#[derive(Debug, Clone, PartialEq)]
pub struct GraphEvent {
    pub(crate) kind: GraphEventType,
    pub(crate) u: u64,
    pub(crate) v: u64,
    pub(crate) ew: f64,
}

impl GraphEvent {
    pub fn new(kind: GraphEventType, u: u64, v: u64, ew: Option<f64>) -> Self {
        Self {
            kind,
            u,
            v,
            ew: ew.unwrap_or(1.),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum GraphEventType {
    NodeAddition = 0,
    NodeRemoval = 1,
    NodeRestoration = 2,
    EdgeAddition = 3,
    EdgeRemoval = 4,
    EdgeWeightUpdate = 5,
    EdgeWeightIncrement = 6,
    TimeStep = 7,
}

impl From<u8> for GraphEventType {
    fn from(value: u8) -> Self {
        match value {
            0 => Self::NodeAddition,
            1 => Self::NodeRemoval,
            2 => Self::NodeRestoration,
            3 => Self::EdgeAddition,
            4 => Self::EdgeRemoval,
            5 => Self::EdgeWeightUpdate,
            6 => Self::EdgeWeightIncrement,
            7 => Self::TimeStep,
            _ => panic!(),
        }
    }
}

pub trait DynAlgorithm {
    fn update(&mut self, e: GraphEvent);
    fn update_batch(&mut self, es: &[GraphEvent]);
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum EdgeDirection {
    InEdges = 0,
    OutEdges = 1,
}
