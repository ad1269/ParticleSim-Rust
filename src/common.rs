use std::time::{SystemTime, UNIX_EPOCH};
use std::io::{BufWriter, Write};
use std::fs::File;

#[macro_export]
macro_rules! max {
    ($x: expr) => ($x);
    ($x: expr, $($z: expr),+) => {{
        let y = max!($($z),*);
        if $x > y {
            $x
        } else {
            y
        }
    }}
}

#[macro_export]
macro_rules! min {
    ($x: expr) => ($x);
    ($x: expr, $($z: expr),+) => {{
        let y = min!($($z),*);
        if $x < y {
            $x
        } else {
            y
        }
    }}
}

// Constants
const DENSITY: f64 = 0.0005;
const MASS: f64 = 0.01;
pub const CUTOFF: f64 = 0.01;
const MIN_R: f64 = (CUTOFF/100.0);
const DT: f64 = 0.0005;

pub const NSTEPS: i32 = 1000;
pub const SAVEFREQ: i32 = 10;

// Global variables
// static mut SIZE: f64 = 0.0;
// static mut FIRST: bool = true;

#[derive(Clone, Copy, Debug)]
pub struct Particle {
	pub x: f64,
	pub y: f64,
	pub vx: f64,
	pub vy: f64,
	pub ax: f64,
	pub ay: f64
}

pub fn get_num_bins(size: f64) -> usize {
	unsafe {
		let bins: usize = (size / CUTOFF).ceil() as usize;
		return bins;
	}
}

pub fn get_time() -> f64 {
	let current_time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap();
	return (current_time.as_secs() as f64) + (current_time.subsec_millis() as f64)/1000.0;
}

pub fn get_size(n: i32) -> f64 {
	return ((n as f64) * DENSITY).sqrt();
}

pub fn init_particles(n: i32, particles: &mut [Particle], size: f64) {
	let sx: i32 = (n as f64).sqrt().ceil() as i32;
	let sy: i32 = (sx + n - 1) / sx;

	let mut shuffle = vec![0; n as usize];
	for i in 0..n {
		shuffle[i as usize] = i;
	}

	for i in 0..n {
		// Make sure particles are not spacially sorted
		let j: usize = rand::random::<usize>() % ((n - i) as usize);
		let k: i32 = shuffle[j];
		shuffle[j] = shuffle[(n - i - 1) as usize];

		// Distributes particles evenly to ensure proper spacing
		// Unsafe because of reference to SIZE
		particles[i as usize].x = size * ((1 + (k % sx)) as f64) / ((1 + sx) as f64);
		particles[i as usize].y = size * ((1 + (k % sx)) as f64) / ((1 + sy) as f64);

		// Assign random velocities within a bound
		particles[i as usize].vx = rand::random::<f64>()*2.0 - 1.0;
		particles[i as usize].vy = rand::random::<f64>()*2.0 - 1.0;
	}
}

pub fn apply_force(particle: &mut Particle, neighbor: &Particle, dmin: &mut f64, davg: &mut f64, navg: &mut i32) {
	let dx = particle.x - neighbor.x;
	let dy = particle.y - neighbor.y;
	let mut r2 = dx*dx + dy*dy;

	if r2 > CUTOFF*CUTOFF {
		return;
	}

	if r2 != 0.0 {
		if r2 / (CUTOFF * CUTOFF) < (*dmin) * (*dmin) {
			*dmin = r2.sqrt() / CUTOFF;
			*davg += r2.sqrt() / CUTOFF;
			*navg += 1;
		}
	}

	r2 = r2.max(MIN_R*MIN_R);
	let r = r2.sqrt();

	let coeff = ( 1.0 - CUTOFF / r ) / r2 / MASS;
	particle.ax += coeff * dx;
	particle.ay += coeff * dy;
}

pub fn move_particle(p: &mut Particle, size: f64) {
	p.vx += p.ax * DT;
	p.vy += p.ay * DT;
	p.x += p.vx * DT;
	p.y += p.vy * DT;

	while p.x < 0.0 || p.x > size {
		if p.x < 0.0 {
			p.x = -p.x
		} else {
			p.x = 2.0*size - p.x
		}
		p.vx = -p.vx
	}

	while p.y < 0.0 || p.y > size {
		if p.y < 0.0 {
			p.y = -p.y
		} else {
			p.y = 2.0*size - p.y
		}
		p.vy = -p.vy
	}
}

pub fn save(writer: &mut BufWriter<File>, n: i32, p: &[Particle], size: f64, is_first: bool) {
	if is_first {
		writeln!(writer, "{} {}", n, size).unwrap();
	}

	for i in 0..n {
		writeln!(writer, "{} {}", p[i as usize].x, p[i as usize].y).unwrap();
	}
}