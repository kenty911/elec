fn main() {
    let file = match netcdf::open("out1.nc") {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Failed to open file: {}", e);
            return;
        }
    };

    println!("Opened file");
    
    // ファイル内の変数のリストを取得
    for var in file.variables() {
        println!("Variable name: {}", var.name());
    }
}
