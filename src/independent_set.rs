use cxx::UniquePtr;

use crate::bridge::{self, *};

pub struct Luby {
    inner: UniquePtr<bridge::Luby>,
}

impl Luby {
    pub fn new() -> Self {
        Self { inner: NewLuby() }
    }
    pub fn run(&mut self, g: &crate::Graph) -> Vec<bool> {
        let mut ret = vec![];
        LubyRun(self.inner.pin_mut(), g, &mut ret);
        ret
    }
    pub fn is_independent_set(&self, set: &[bool], g: &crate::Graph) -> bool {
        LubyIsIndependentSet(&self.inner, set, g)
    }
}
