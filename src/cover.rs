use std::ops::Deref;

use cxx::UniquePtr;
use miette::{bail, IntoDiagnostic, Result};

use crate::bridge::{self, *};

pub struct Cover {}
