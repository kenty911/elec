fn main() {
    println!("cargo:rustc-link-lib=netcdf");
    println!("cargo:rustc-link-search=native=/opt/local/netcdf/lib");
}
