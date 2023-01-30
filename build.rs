use std::env;
use std::path::PathBuf;

fn main() {
    let mut b = autocxx_build::Builder::new("src/main.rs", &[
        "networkit/include",
        "networkit/extlibs/tlx",
        "bridge/",
        "/opt/homebrew/opt/libomp/include"
    ]).build().unwrap();
    // This assumes all your C++ bindings are in main.rs
    b.flag_if_supported("-std=c++17")
        .compile("networkit-rs"); // arbitrary library name, pick anything
    println!("cargo:rerun-if-changed=src/main.rs");

    // Add instructions to link to any C++ libraries you need.

    build_networkit();
}

fn build_networkit() {
    // let base = env::var_os("OUT_DIR").unwrap();
    // let mut out_lib_path = PathBuf::from(base);
    // out_lib_path.push("libnetworkit.a");
    // if out_lib_path.exists() {
    //     return;
    // }

    let mut config = cmake::Config::new("networkit");
    config.define("NETWORKIT_STATIC", "ON");
    let mut res = config.build();
    println!("{}", res.display());
    println!("cargo:rustc-link-search=native={}", res.display());
    res.push("lib");
    println!("cargo:rustc-link-search=native={}", res.display());
    println!("cargo:rustc-link-lib=static=networkit");
}