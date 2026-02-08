use mesh::generate_mesh;
use problem_loader::load_problem_from_file;
use pde::*;
use util::gauss_pp::gauss_pp;
mod plotter;
mod mesh;
mod pde;
mod problem_loader;

use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
struct Cli {
    /// The path to the file to read
    path: PathBuf,
}

fn main() {
    println!("Hello, world!");
    // TODO Add problem stmt
    let  args = Cli::parse();

    // Abstracted loading
    let problem = load_problem_from_file(&args.path);
    let soln_function = |x:f64| (x*((1 as f64)-x))/2 as f64;
    
    println!("Loaded Parameters: p={}, q={}, f={}", 
             problem.p_val, problem.q_val, problem.f_val);
    println!("Begin SIPDG Process");

    // ------------------- Generate Mesh --------------------
    // Domain to find soln for
    let domain_a: f64 = 0 as f64; 
    let domain_b: f64 = 1.0;

    // number of elements
    let num_elements: usize = 30;

    // Generate the mesh
    let (h_elem, x_dof) = generate_mesh(domain_a, domain_b, num_elements);
    
    println!("Elements array \n {:?}", h_elem);
    println!("x dof is \n {:?}", x_dof);

    let n_dof: usize = x_dof.len();
    let mut assembler = pde::SipdgAssembler::new(h_elem, x_dof, 80.0);

    // -------------------- Assemble system matrix ---------
    assembler.assemble_volume(&problem);
    assembler.assemble_interfaces(&problem);

    // Define the specific conditions for this run
    let left_bc = DirichletBC { value: 0.0 };  // u(0) = 0
    let right_bc = DirichletBC { value: 0.0 }; // u(1) = 0

    assembler.apply_boundaries(&problem, &left_bc, &right_bc);
                                         
    // solve system
    let u = gauss_pp(assembler.a, assembler.rhs);

    let mut l2_err_sq: f64 = 0.0;
    let npts: usize = 10;

    for k in 0..num_elements {
        let idx0 = 2 * k;
        let idx1 = 2 * k + 1;

        let x_l = assembler.x_dof[idx0];
        let x_r = assembler.x_dof[idx1];
        let h = assembler.h_elem[k];

        let ul = u[idx0];
        let ur = u[idx1];

        let step = (x_r - x_l) / ((npts - 1) as f64);

        // trapezoid integration on-the-fly
        let mut prev_x = x_l;
        let mut prev_err2 = {
            let phi1 = (x_r - prev_x) / h;
            let phi2 = (prev_x - x_l) / h;
            let uh = ul * phi1 + ur * phi2;
            let uex = soln_function(prev_x);
            let e = uh - uex;
            e * e
        };

        for i in 1..npts {
            let x = x_l + step * (i as f64);

            let phi1 = (x_r - x) / h;
            let phi2 = (x - x_l) / h;

            let uh = ul * phi1 + ur * phi2;
            let uex = soln_function(x);
            let e = uh - uex;
            let err2 = e * e;

            let dx = x - prev_x;
            l2_err_sq += 0.5 * dx * (prev_err2 + err2);

            prev_x = x;
            prev_err2 = err2;
        }
    }

    let l2_error = l2_err_sq.sqrt();
    println!("Total DoFs: {}", n_dof);
    println!("L2 Error: {:.6e}", l2_error);

    if let Err(e) = plotter::plot_results(&assembler.x_dof, &u, soln_function, n_dof) {
        eprintln!("Failed to create plot: {}", e);
    }
}
