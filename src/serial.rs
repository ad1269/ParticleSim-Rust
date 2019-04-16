#[path = "common.rs"]
#[macro_use]
mod common;

use std::fs::File;
use std::io::BufWriter;
use getopts::Options;
use std::env;

use common::{Particle, CUTOFF, NSTEPS, SAVEFREQ};

pub fn simulate_main() {
	let (mut navg, mut nabsavg): (i32, i32) = (0, 0);
	let (mut absmin, mut absavg, mut davg, mut dmin): (f64, f64, f64, f64) = (1.0, 1.0, 0.0, 0.0);

	// Set and parse command line options
	let args: Vec<String> = env::args().collect();
	let mut opts = Options::new();
	opts.optopt("n", "", "set number of particles", "NUM");
    opts.optopt("o", "", "set output file name", "NAME");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!(f.to_string()) }
    };

    // Read and set command line options
	let n: i32 = if matches.opt_present("n") {
		matches.opt_str("n").unwrap().parse::<i32>().unwrap()
	} else {
		500
	};

	let savename: String = if matches.opt_present("o") {
		matches.opt_str("o").unwrap()
	} else {
		"out".to_string()
	};

	// Open output file and writer
	let write_file = File::create(savename).unwrap();
    let mut writer = BufWriter::new(&write_file);

    // Set up binning
    let bins: usize = common::get_num_bins();
    let mut pointers: Vec<Vec<Particle>> = Vec::new();
    for i in &mut pointers {
    	*i = Vec::new();
    }

    // Populate bins
    {
    	// Initialize particles
	    let mut particles: Box<[Particle]> = vec![Particle{x: 0., y: 0., vx: 0., vy: 0., ax: 0., ay: 0.}; n as usize].into_boxed_slice();
	    common::set_size(n);
	    common::init_particles(n, &mut particles);

	    // Categorize into bins
	    for i in 0..n {
	    	let bin_i: usize = min!((particles[i as usize].x / CUTOFF) as usize, bins - 1);
	    	let bin_j: usize = min!((particles[i as usize].y / CUTOFF) as usize, bins - 1);

	        pointers[bin_i*bins + bin_j].push(particles[i as usize]);
	    }
    }
    

    let mut simulation_time = common::get_time();
    for step in 0..NSTEPS {
    	navg = 0;
    	davg = 0.;
    	dmin = 1.;

    	// Compute forces
    	for bin_i in 0..bins {
    		for bin_j in 0..bins {
    			for i in 0..pointers[bin_i*bins + bin_j].len() {
		    		pointers[bin_i*bins + bin_j][i].ax = 0.;
		    		pointers[bin_i*bins + bin_j][i].ay = 0.;

		    		for other_bin_i in (max!(0, (bin_i as i32) - 1) as usize)..=min!(bins - 1, bin_i + 1) {
		    			for other_bin_j in (max!(0, (bin_j as i32) - 1) as usize)..=min!(bins - 1, bin_j + 1) {
		    				for k in 0..pointers[other_bin_i*bins + other_bin_j].len() {
		    					common::apply_force(&mut pointers[bin_i*bins + bin_j][i], &pointers[other_bin_i*bins + other_bin_j][k], &mut dmin, &mut davg, &mut navg);
		    				}
		    			}
		    		}
		    	}
	    	}
    	}

    }
}