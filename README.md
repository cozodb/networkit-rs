# NetworKit-rs

[![Crates.io](https://img.shields.io/crates/v/networkit-rs)](https://crates.io/crates/networkit-rs)
[![docs.rs](https://img.shields.io/docsrs/networkit-rs)](https://docs.rs/networkit-rs)


This is a rust wrapper for [NetworKit](https://networkit.github.io/), a toolkit for large-scale network analysis.

The API is similar to the python API of the upstream package.

This is still a work in progress, hence incomplete. If you need some function that is not available, please open an issue.

## Safety issues

The upstream library does not check for out-of-bound array access most of the time and leads to undefined behaviour when out-of-bound indices are given as arguments. This occurs for both the C++ code and the official Python wrapper. For the current rust wrapper, it is simply infeasible to rectify this for all usage (though we did try to patch it up for some usage). So you **must** ensure that indices are valid in your code, even when calling functions that are not marked unsafe.

In addition, the upstream C++ library may throw exceptions, which are not always handled on the Rust side yet (since there are simply too many exceptions thrown, and currently every one of them must be handled manually). So for some of them, the Rust code will simply panic. If you encounter such problems, open an issue.

Due to the above problems, it is not recommended to use this library for user-facing programs. Using it internally in well-tested environments where every input is validated is probably OK.