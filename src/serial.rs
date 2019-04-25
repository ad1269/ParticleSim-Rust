#[path = "common.rs"]
#[macro_use]
mod common;

use std::fs::File;
use std::io::BufWriter;
use getopts::Options;
use std::env;

use common::{Particle, CUTOFF, NSTEPS, SAVEFREQ};

pub fn simulate_main() {
	let mut nabsavg: i32 = 0;
	let mut navg: i32;
	let (mut absmin, mut absavg): (f64, f64) = (1.0, 1.0);
 	let (mut davg, mut dmin): (f64, f64);

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
    let mut writer: BufWriter<File> = BufWriter::new(write_file);

    // Set up binning
    let size: f64 = common::get_size(n);
    let bins: usize = common::get_num_bins(size);
    let mut pointers: Vec<Vec<usize>> = vec![Vec::new(); bins * bins];
    
	  // Initialize particles
    let mut particles = vec![Particle{x: 0., y: 0., vx: 0., vy: 0., ax: 0., ay: 0.}; n as usize];
    common::init_particles(n, &mut particles, size);

    // Categorize into bins
    for i in 0..n {
    	let bin_i: usize = min!((particles[i as usize].x / CUTOFF) as usize, bins - 1);
    	let bin_j: usize = min!((particles[i as usize].y / CUTOFF) as usize, bins - 1);

      pointers[bin_i*bins + bin_j].push(i as usize);
    }

    let mut is_first_save: bool = true;
    let mut simulation_time = common::get_time();
    for step in 0..NSTEPS {
    	let mut navg = 0;
    	let mut davg = 0.;
    	let mut dmin = 1.;
        //
        //  compute forces
        //
        for i in 0..particles.len() {
            particles[i as usize].ax = 0.;
            particles[i as usize].ay = 0.;

	    	let bin_i: usize = min!((particles[i as usize].x / CUTOFF) as usize, bins - 1);
	    	let bin_j: usize = min!((particles[i as usize].y / CUTOFF) as usize, bins - 1);
	    	// let index0: usize = pointers[bin_i*bins + bin_j][i];

    		for other_bin_i in (max!(0, (bin_i as i32) - 1) as usize)..=(min!(bins - 1, bin_i + 1) as usize) {
    			for other_bin_j in (max!(0, (bin_j as i32) - 1) as usize)..=(min!(bins - 1, bin_j + 1) as usize) {
               		for k in 0..pointers[other_bin_i*bins + other_bin_j].len() {
               			let index1: usize = pointers[other_bin_i*bins + other_bin_j][k];
               			let particle1: Particle = particles[index1];
	    				common::apply_force(&mut particles[i], &particle1, &mut dmin, &mut davg, &mut navg);
	    			}
                }
            }
        }
 
        //
        //  move particles
        //
        for i in 0..particles.len() {

	    	let old_bin_i: usize = min!((particles[i as usize].x / CUTOFF) as usize, bins - 1);
	    	let old_bin_j: usize = min!((particles[i as usize].y / CUTOFF) as usize, bins - 1);

            common::move_particle( &mut particles[i], size ); 

	    	let new_bin_i: usize = min!((particles[i as usize].x / CUTOFF) as usize, bins - 1);
	    	let new_bin_j: usize = min!((particles[i as usize].y / CUTOFF) as usize, bins - 1);

            // Has been moved into a new bin
            if old_bin_j != new_bin_j || old_bin_i != new_bin_i {

            	let old_bin = &mut pointers[old_bin_i*bins + old_bin_j];
				let index = old_bin.iter().position(|x| *x == i).unwrap();
				old_bin.remove(index);

				let new_bin = &mut pointers[new_bin_i*bins + new_bin_j];
				new_bin.push(i);
            }
        }
        
          //
          // Computing statistical data
          //
          if navg != 0 {
            absavg +=  davg/(navg as f64);
            nabsavg += 1;
          }
          if dmin < absmin {
          	absmin = dmin;
          }
        
          //
          //  save if necessary
          //
          if (step%SAVEFREQ) == 0 {
              common::save( &mut writer, n, &particles, size, is_first_save );
              is_first_save = false;
          }
    }
    simulation_time = common::get_time() - simulation_time;
    print!( "n = {}, simulation time = {} seconds", n, simulation_time);


    if nabsavg != 0 {
    	absavg /= nabsavg as f64;
    }

    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    println!(", absmin = {}, absavg = {}", absmin, absavg);
    if absmin < 0.4 {
    	println! ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    }
    if absavg < 0.8 {
    	println! ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
}