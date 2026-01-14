use math::linespace::linespace;
use math::matrix::Matrix;

fn main() {
    println!("Hello, world!");


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
    let A:Matrix = Matrix::new(n_dof, n_dof);     
    let F:Vec<f64> = vec![0.0; n_dof];

    for i in 0..num_elements {
        // element parameters
        let h_k = h_elem[i];
        let idx = [2*i, 2*i+1];
        //let xc = (x_dof[idx[1]] + x_dof[idx[2]]) / 2.0;


        // Evaluate Coeffs at midpoint
        //let p_k = p_func(xc);
        //let q_k = q_func(xc);
        //let f_k = f_func(xc);

        //  % Gradients of basis functions on reference element [-1, 1] mapped to h
        // dphi/dx = +/- 1/h. 
        // Local Stiffness Matrix (Int p u' v')
        // integral (-1/h)(-1/h) dx = 1/h^2 * h = 1/h 
        //let k_e: Matrix = (p_k / h_k) * Matrix::new(2, 2);

    }

}
