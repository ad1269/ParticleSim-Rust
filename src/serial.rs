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
    let mut writer: BufWriter<File> = BufWriter::new(write_file);

    // Set up binning
    let size: f64 = common::get_size(n);
    let bins: usize = common::get_num_bins(size);
    let mut pointers: Vec<Vec<usize>> = vec![Vec::new(); bins * bins];

    // Populate bins
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
    	navg = 0;
    	davg = 0.;
    	dmin = 1.;
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

            	let mut old_bin = &mut pointers[old_bin_i*bins + old_bin_j];
				let index = old_bin.iter().position(|x| *x == i).unwrap();
				old_bin.remove(index);

				let mut new_bin = &mut pointers[new_bin_i*bins + new_bin_j];
				new_bin.push(i);

                // std::vector<particle_t*> &old_bin = pointers[old_bin_i*bins + old_bin_j];
                // old_bin.erase(std::find(old_bin.begin(), old_bin.end(), &particles[i]));

                // std::vector<particle_t*> &new_bin = pointers[new_bin_i*bins + new_bin_j];
                // new_bin.push_back(&particles[i]);
            }
        }


        // if( find_option( argc, argv, "-no" ) == -1 )
        // {
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

        // }
    }
   

    // let mut simulation_time = common::get_time();
    // for step in 0..NSTEPS {
    // 	navg = 0;
    // 	davg = 0.;
    // 	dmin = 1.;

    // 	// Compute forces
    // 	for bin_i in 0..bins {
    // 		for bin_j in 0..bins {
    // 			for i in 0..pointers[bin_i*bins + bin_j].len() {
		  //   		pointers[bin_i*bins + bin_j][i].ax = 0.;
		  //   		pointers[bin_i*bins + bin_j][i].ay = 0.;

		  //   		for other_bin_i in (max!(0, (bin_i as i32) - 1) as usize)..=min!(bins - 1, bin_i + 1) {
		  //   			for other_bin_j in (max!(0, (bin_j as i32) - 1) as usize)..=min!(bins - 1, bin_j + 1) {
		  //   				for k in 0..pointers[other_bin_i*bins + other_bin_j].len() {
		  //   					common::apply_force(&mut pointers[bin_i*bins + bin_j][i], &pointers[other_bin_i*bins + other_bin_j][k], &mut dmin, &mut davg, &mut navg);
		  //   				}
		  //   			}
		  //   		}
		  //   	}
	   //  	}
    // 	}

    // }
}