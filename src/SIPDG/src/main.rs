use util::linespace::linespace;
use util::matrix::Matrix;
use util::gauss_pp::gauss_pp;

use clap::Parser;
use std::path::PathBuf;
use std::fs;

#[derive(Parser)]
#[command(name = "Reader")]
#[command(about = "Reads a specific file")]
struct Cli {
    /// The path to the file to read
    path: PathBuf,

    /// Optional: A flag to print the word count
    #[arg(short, long)]
    count: bool,
}

fn main() {
    println!("Hello, world!");
    // TODO Add problem stmt
    let cli = Cli::parse();

    // Try to read the file content
    match fs::read_to_string(&cli.path) {
        Ok(content) => {
            if cli.count {
                println!("File contains {} words.", content.split_whitespace().count());
            } else {
                println!("--- File Content ---\n{}", content);
            }
        }
        Err(e) => {
            eprintln!("Error reading {:?}: {}", cli.path, e);
        }
    }


    println!("Begin SIPDG Process");

    // Problem: -(p(x)u')' + q(x)u = f(x) on [0,1]
    // Note these could have to be lambdas in the future
    let p_func = |x: f64| 1.0;
    let q_func = |x: f64| 0.0;
    let f_func = |x: f64| 1.0;
    let soln_function = |x:f64| (x*((1 as f64)-x))/2 as f64;
    
    // Penalty parameter for stability
    let penalty_param: u32 = 10;

    // ------------------- Generate Mesh --------------------

    // Domain to find soln for
    let domain_a: f64 = 0 as f64; 
    let domain_b: f64 = 1.0;

    // number of elements
    let num_elements: usize = 20;

    // Generate evenly spaced points across the domain
    let x_interface: Vec<f64> = linespace(domain_a, domain_b, num_elements + 1);
    println!("Evenly spaced points are \n {:?}", x_interface);

    // With DG we do not share nodes 
    let n_dof: usize = 2 * num_elements;
    let mut x_dof: Vec<f64> = vec![0.0; n_dof];
    let mut h_elem: Vec<f64> = vec![0.0; num_elements];

    // Fill node coordinates and element sizes
    for i in 0..num_elements {
        let x_l: f64 = x_interface[i];
        let x_r: f64 = x_interface[i + 1];
        h_elem[i] = x_r - x_l;

        // map global DoF indices for this element
        x_dof[2*i] = x_l;
        x_dof[2*i+1] = x_r;
    }
    println!("Elements array \n {:?}", h_elem);
    println!("x dof is \n {:?}", x_dof);

    // -------------------- Assemble system matrix ---------
    let mut A:Matrix<f64> = Matrix::<f64>::new(n_dof, n_dof, 0.0);     
    let mut F:Matrix<f64> = Matrix::<f64>::new(1, n_dof, 0.0);

    for i in 0..num_elements {
        // element parameters
        let h_k = h_elem[i];
        let idx = [2*i, 2*i+1];
        let xc = (x_dof[idx[0]] + x_dof[idx[1]]) / 2.0;


        // Evaluate Coeffs at midpoint
        let p_k = p_func(xc);
        let q_k = q_func(xc);
        let f_k = f_func(xc);

        //  % Gradients of basis functions on reference element [-1, 1] mapped to h
        // dphi/dx = +/- 1/h. 
        // Local Stiffness Matrix (Int p u' v')
        // integral (-1/h)(-1/h) dx = 1/h^2 * h = 1/h 
        let mut multiplication_matrix = Matrix::<f64>::new(2, 2, -1.0);
        multiplication_matrix[(0, 0)] = 1.0;
        multiplication_matrix[(1, 1)] = 1.0;
        let k_e: Matrix<f64> = &multiplication_matrix * &(p_k / h_k) ;

        // Local Mass Matrix (Int q u v)
        // Standard linear mass matrix [2 1; 1 2] * h/6
        multiplication_matrix[(0, 0)] = 2.0;
        multiplication_matrix[(0, 1)] = 1.0;
        multiplication_matrix[(1, 0)] = 1.0;
        multiplication_matrix[(1, 1)] = 2.0;
        let m_e: Matrix<f64> = &multiplication_matrix * &(q_k * h_k / 6.0);


        // local load vector
        let multi_matrix = Matrix::<f64>::new(1, 2, 1.0);
        let f_e = &multi_matrix * &(f_k * h_k / 2.0);

        // assemble volume terms
        for a in 0..2 {
            let I = idx[a];
            F[(0,I)] += f_e[(0,a)];

            for b in 0..2 {
                let J = idx[b];
                A[(I, J)] += k_e[(a, b)] + m_e[(a, b)];
            }
        }

    }

    // interface integral
    for i in 0..num_elements-1 {
        // interface is between R_i and L_{i+1}
        let indx_l: usize = 2 * i + 1;     // right node of element i
        let indx_r: usize = 2 * i + 2;     // left node of element i+1

        let x_val = x_dof[indx_l];         // interface x
        let p_val = p_func(x_val);

        let h_avg = 0.5 * (h_elem[i] + h_elem[i + 1]);
        let penalty = penalty_param as f64 * (p_val / h_avg);

        let grad_phi_l = 1.0 / h_elem[i];       // derivative on element i at its right end
        let grad_phi_r = -1.0 / h_elem[i + 1];  // derivative on element i+1 at its left end
                                                
        A[(indx_l, indx_l)] = A[(indx_l, indx_l)] + penalty;
        A[(indx_l, indx_r)] = A[(indx_l, indx_r)] - penalty;
        A[(indx_r, indx_l)] = A[(indx_r, indx_l)] - penalty;
        A[(indx_r, indx_r)] = A[(indx_r, indx_r)] + penalty;

        A[(indx_l, indx_l)] = A[(indx_l, indx_l)] - 0.5 * p_val * (1.0/h_elem[i]);
        A[(indx_l, indx_r)] = A[(indx_l, indx_r)] - 0.5 * p_val * (-1.0/h_elem[i+1]);

        A[(indx_l, indx_l)] = A[(indx_l, indx_l)] - 0.5 * p_val * (1.0/h_elem[i]);
        A[(indx_r, indx_l)] = A[(indx_r, indx_l)] + 0.5 * p_val * (1.0/h_elem[i]);

        
        A[(indx_r, indx_l)] = A[(indx_r, indx_l)] + 0.5 * p_val * (1.0/h_elem[i]);
        A[(indx_r, indx_r)] = A[(indx_r, indx_r)] + 0.5 * p_val * (-1.0/h_elem[i+1]);

        A[(indx_l, indx_r)] = A[(indx_l, indx_r)] - 0.5 * p_val * grad_phi_r + 0.5 * p_val * grad_phi_l; 
        A[(indx_r, indx_l)] = A[(indx_r, indx_l)] + 0.5 * p_val * grad_phi_l - 0.5 * p_val * grad_phi_r; 

    }

    let bnd_nodes: [usize; 2] = [0, n_dof - 1];
    let normals: [f64; 2] = [-1.0, 1.0];

    for side in 0..2 {
        let idx = bnd_nodes[side]; // global dof index
        let n = normals[side];     // outward normal

        // Pick boundary element size and basis gradient at the boundary node
        let (h_bnd, grad_phi) = if side == 0 {
            // Left boundary: first element, left node basis slope is -1/h
            let h = h_elem[0];
            (h, -1.0 / h)
        } else {
            // Right boundary: last element, right node basis slope is +1/h
            let h = h_elem[num_elements - 1];
            (h, 1.0 / h)
        };

        let p_val = p_func(x_dof[idx]);

        // Penalty term: sigma_0 * p / h
        let pen = (penalty_param as f64) * (p_val / h_bnd);
        A[(idx, idx)] += pen;

        // Consistency + symmetry terms:
        A[(idx, idx)] -= 2.0 * p_val * grad_phi * n;

    }

    // solve system
    let u = gauss_pp(A, F);

    let mut l2_err_sq: f64 = 0.0;
    let npts: usize = 10;

    for k in 0..num_elements {
        let idx0 = 2 * k;
        let idx1 = 2 * k + 1;

        let x_l = x_dof[idx0];
        let x_r = x_dof[idx1];
        let h = h_elem[k];

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

}
