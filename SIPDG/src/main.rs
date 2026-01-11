use math::linespace::linespace;

fn main() {
    println!("Hello, world!");


    println!("Begin SIPDG Process");

    // Problem: -(p(x)u')' + q(x)u = f(x) on [0,1]
    // Note these could have to be lambdas in the future
    let p_coeff: u32 = 1;
    let q_coeff: u32 = 0;
    let f_coeff: u32 = 1;
    let soln_function = |x:f64| (x*((1 as f64)-x))/2 as f64;
    
    println!("Our function looks like -({}u')' + {}u = {}", p_coeff, q_coeff, f_coeff);


    // Penatly parameter for stability
    let penalty_param: u32 = 10;

    // ------------------- Generate Mesh --------------------

    // Domain to find soln for
    let domain_a: f64 = 0 as f64; 
    let domain_b: f64 = 1 as f64;

    // number of elements
    let num_elements: usize = 20;

    // Generate evenly spaced points across the domain
    let x_interface: Vec<f64> = linespace(domain_a, domain_b, num_elements);
    println!("Evenly spaced points are \n {:?}", x_interface);

    // With DG we do not share nodes 
    let n_dof: usize = 2 * num_elements;
    let mut x_dof: Vec<f64> = vec![0 as f64; n_dof];
    let mut h_elem: Vec<f64> = vec![0 as f64; num_elements];

    // Fill node coordinates and element sizes

    
}
