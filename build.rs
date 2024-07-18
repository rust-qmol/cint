use cmake::Config;
use std::process::Command;

fn main() {
    Command::new("git")
        .args(&["submodule", "init"])
        .status()
        .expect("git submodule init failed");
    Command::new("git")
        .args(&["submodule", "update"])
        .status()
        .expect("git submodule update faild");

    let dst = Config::new("libcint").build();

    println!("cargo:rustc-link-search=native={}/lib", dst.display());
    println!("cargo:rustc-link-lib=dylib=cint");

    let bindings = bindgen::Builder::default()
        .clang_args(&[format!("-I{}/include/", dst.display())])
        .header_contents("const_cint.h", "
        #define HAVE_DEFINED_CINTINTEGRALFUNCTION
        
        typedef void CINTOptimizerFunction(CINTOpt const **opt,
            FINT const atm[], FINT natm, FINT const bas[], FINT nbas, double const env[]);
        
        typedef CACHE_SIZE_T CINTIntegralFunction(double out[], FINT const dims[], FINT const shls[],
            FINT const atm[], FINT natm, FINT const bas[], FINT nbas, double const env[],
            CINTOpt const *opt, double cache[]);
        
        ")
        .header(format!("{}/include/cint.h", dst.display()))
        .header(format!("{}/include/cint_funcs.h", dst.display()))
        .allowlist_function("int.*")
        .allowlist_function(".*_optimizer")
        .allowlist_function("cint2e.*")
        .allowlist_function("CINT.*")
        .allowlist_type("CINT.*")
        .allowlist_var(".*")
        .allowlist_type("CINTOptimizerFunction")
        .allowlist_type("CINTIntegralFunction")
        .blocklist_item("_.*")
        .dynamic_link_require_all(true)
        .wrap_unsafe_ops(true)
        .generate_comments(true)
        .layout_tests(false)
        .generate()
        .expect("Unable to generate bindings");

    let out_path = std::env::current_dir().unwrap();
    bindings
        .write_to_file(out_path.join("src/cint/libcint.rs"))
        .expect("Couldn't write bindings!");
}
