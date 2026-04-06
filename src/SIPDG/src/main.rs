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
    /// The basis function degree (1 for linear, 2 for quadratic)
    #[arg(short, long, default_value_t = 1)]
    degree: usize,
}

fn main() {
    println!("Hello, world!");
    // TODO Add problem stmt
    let args = Cli::parse();

    let order = match args.degree {
        1 => BasisOrder::Linear,
        2 => BasisOrder::Quadratic,
        _ => panic!("Unsupported degree: {}. Use 1 or 2.", args.degree),
    };

    // Abstracted loading
    let problem = load_problem_from_file(&args.path);
    let soln_function = |x: f64| (x * (1.0 - x)) / 2.0;
    
    println!("Loaded Parameters: a={}, q={}, f={}", 
             problem.a_val, problem.q_val, problem.f_val);
    println!("Basis Order: {:?}", order);
    println!("Begin SIPDG Process");

    // ------------------- Generate Mesh --------------------
    // Domain to find soln for
    let domain_a: f64 = 0.0; 
    let domain_b: f64 = 1.0;

    // number of elements
    let num_elements: usize = 3000;

    // Generate the mesh
    let (h_elem, x_dof) = generate_mesh(domain_a, domain_b, num_elements, order);
    
    // println!("Elements array \n {:?}", h_elem);
    // println!("x dof is \n {:?}", x_dof);

    let n_dof: usize = x_dof.len();
    let mut assembler = pde::SipdgAssembler::new(h_elem, x_dof, 20.0, order);

    // -------------------- Assemble system matrix ---------
    assembler.assemble_volume(&problem);
    assembler.assemble_interfaces(&problem);

    let (mut a, mut rhs) = assembler.assemble_to_global();

    // Define the specific conditions for this run
    let left_bc = DirichletBC { value: 0.0 };  // u(0) = 0
    let right_bc = DirichletBC { value: 0.0 }; // u(1) = 0

    assembler.apply_boundaries(&problem, &left_bc, &right_bc, &mut a, &mut rhs);
                                         
    // solve system
    let p = gauss_pp(a, rhs);

    let mut l2_err_sq: f64 = 0.0;
    let quad_pts = util::quadrature::get_gauss_legendre_3pts();

    for elem in &assembler.elements {
        let det_j = elem.jacobian();
        for pt in &quad_pts {
            let x = elem.map_to_physical(pt.xi);
            let mut u_h = 0.0;
            for i in 0..order.num_nodes() {
                u_h += p[elem.nodes[i]] * elem.phi(i, pt.xi);
            }
            let u_ex = soln_function(x);
            let err = u_h - u_ex;
            l2_err_sq += pt.weight * err * err * det_j;
        }
    }

    let l2_error = l2_err_sq.sqrt();
    println!("Total DoFs: {}", n_dof);
    println!("L2 Error: {:.6e}", l2_error);

    if let Err(e) = plotter::plot_results(&assembler.x_dof, &p, soln_function, num_elements, order) {
        eprintln!("Failed to create plot: {}", e);
    }
}
