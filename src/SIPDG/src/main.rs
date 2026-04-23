use mesh::generate_mesh;
use problem_loader::load_problem_from_file;
use pde::*;
use util::cg::cg;
use rayon::prelude::*;
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
    
    // Exact solution for the problem: - (a p')' + q p = f 
    let soln_function: Box<dyn Fn(f64) -> f64 + Send + Sync> = match problem.problem_type {
        ProblemType::Sine => {
            println!("Using SL Sine Exact Solution: sin(pi * x)");
            Box::new(|x: f64| (std::f64::consts::PI * x).sin())
        }
        ProblemType::Cosine => {
            println!("Using SL Cosine Exact Solution: cos(pi * x)");
            Box::new(|x: f64| (std::f64::consts::PI * x).cos())
        }
        ProblemType::Constant => {
            let a_val = problem.a_val;
            let q_val = problem.q_val;
            let f_val = problem.f_val;

            if q_val > 0.0 {
                println!("Using General SL Exact Solution (q > 0)");
                let sqrt_qa = (q_val / a_val).sqrt();
                Box::new(move |x: f64| {
                    (f_val / q_val) * (1.0 - (sqrt_qa * (x - 0.5)).cosh() / (0.5 * sqrt_qa).cosh())
                })
            } else {
                println!("Using Poisson Exact Solution (q = 0)");
                Box::new(move |x: f64| (f_val / a_val) * (x * (1.0 - x)) / 2.0)
            }
        }
    };
    
    println!("Loaded Parameters: a={}, q={}, f={}, num_elements={}, sigma_0={}", 
             problem.a_val, problem.q_val, problem.f_val, problem.num_elements, problem.sigma_0);
    println!("Basis Order: {:?}", order);
    println!("Begin SIPDG Process");

    // ------------------- Generate Mesh --------------------
    // Domain to find soln for
    let domain_a: f64 = 0.0; 
    let domain_b: f64 = 1.0;

    // number of elements
    let num_elements = problem.num_elements;

    // Generate the mesh
    let (h_elem, x_dof) = generate_mesh(domain_a, domain_b, num_elements, order);
    
    // println!("Elements array \n {:?}", h_elem);
    // println!("x dof is \n {:?}", x_dof);

    let n_dof: usize = x_dof.len();
    let mut assembler = pde::SipdgAssembler::new(h_elem, x_dof, problem.sigma_0, order);

    // -------------------- Assemble system matrix ---------
    assembler.assemble_volume(&problem);
    assembler.assemble_interfaces(&problem);

    let (mut a, mut rhs) = assembler.assemble_to_global();

    // Define the specific conditions for this run
    let (left_bc, right_bc) = match problem.problem_type {
        ProblemType::Cosine => {
            // p(0) = 1 (Dirichlet), p'(1) = 0 (Neumann)
            (BCType::Dirichlet(DirichletBC { value: 1.0 }), BCType::Neumann(NeumannBC { value: 0.0 }))
        }
        _ => {
            // Default: Homogeneous Dirichlet p(0) = 0, p(1) = 0
            (BCType::Dirichlet(DirichletBC { value: 0.0 }), BCType::Dirichlet(DirichletBC { value: 0.0 }))
        }
    };

    assembler.apply_boundaries(&problem, &left_bc, &right_bc, &mut a, &mut rhs);
                                         
    // solve system (Matrix-Free)
    let op = assembler.matrix_free_op(&problem, &left_bc, &right_bc);
    let p = cg(&op, &rhs, 1e-10, 10000);

    let quad_pts = util::quadrature::get_gauss_legendre_3pts();
    let n_nodes = order.num_nodes();

    let l2_err_sq: f64 = assembler.elements.par_iter().map(|elem| {
        let det_j = elem.jacobian();
        let mut local_err_sq = 0.0;
        for pt in &quad_pts {
            let x = elem.map_to_physical(pt.xi);
            let mut u_h = 0.0;
            for i in 0..n_nodes {
                u_h += p[elem.nodes[i]] * elem.phi(i, pt.xi);
            }
            let u_ex = soln_function(x);
            let err = u_h - u_ex;
            local_err_sq += pt.weight * err * err * det_j;
        }
        local_err_sq
    }).sum();

    let l2_error = l2_err_sq.sqrt();
    println!("Total DoFs: {}", n_dof);
    println!("L2 Error: {:.6e}", l2_error);

    if let Err(e) = plotter::plot_results(&assembler.x_dof, &p, soln_function, num_elements, order) {
        eprintln!("Failed to create plot: {}", e);
    }
}
