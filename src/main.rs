mod serial;
mod parallel;

use std::env;
use std::fs::File;
use std::io::BufWriter;
use getopts::Options;

fn main() {
    // Set and parse command line options
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    opts.optopt("p", "", "set parallel mode", "PARALLEL");
    opts.optopt("n", "", "set number of particles", "NUM");
    opts.optopt("o", "", "set output file name", "NAME");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!(f.to_string()) }
    };

    // Read and set command line options
    // Parallel vs serial mode
    let is_parallel: bool = if matches.opt_present("p") {
        matches.opt_str("p").unwrap().parse::<i32>().unwrap() != 0
    } else {
    	panic!("must provide parallel or serial mode");
    };
    println!("parallel mode: {:}", is_parallel);
    // Num particles
    let n: i32 = if matches.opt_present("n") {
        matches.opt_str("n").unwrap().parse::<i32>().unwrap()
    } else {
        500
    };
    // Output file
    let savename: String = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        "out".to_string()
    };

    // Open output file and writer
    let write_file = File::create(savename).unwrap();
    let mut writer: BufWriter<File> = BufWriter::new(write_file);

    if is_parallel {
    	parallel::simulate_main(&mut writer, n);
    } else {
    	serial::simulate_main(&mut writer, n);
    }
}
