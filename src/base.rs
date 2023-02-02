use miette::Result;

pub trait Algorithm {
    fn run(&mut self) -> Result<()>;
    fn has_finished(&self) -> bool;
}
