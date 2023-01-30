use std::env;

fn main() {
    let target = env::var("TARGET").unwrap();

    let mut builder = cxx_build::bridge("src/bridge.rs");
    if target.contains("msvc") {
        builder.flag_if_supported("-std:c++17");
    } else {
        builder.flag("-std=c++17");
    };

    builder
        .files(["bridge/bridge.cpp"])
        .include("networkit/include")
        .include("networkit/extlibs/tlx")
        .include("bridge")
        .include("/opt/homebrew/opt/libomp/include");

    println!("cargo:rerun-if-changed=src/bridge.rs");
    println!("cargo:rerun-if-changed=bridge/bridge.h");

    println!("cargo:rustc-link-search=native=/opt/homebrew/opt/libomp/lib");
    println!("cargo:rustc-link-lib=static=omp");
    builder.compile("networkit-rs");

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