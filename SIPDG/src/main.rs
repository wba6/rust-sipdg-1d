
/// generates a row vector of a specified number of linearly (evenly) spaced points between two endpoints
fn linespace(start: f64, end: f64, num_elements: usize) -> Vec<f64> {
    assert!(num_elements > 0);
    // Create a vector or num_elements + 1 due to this being inclusive of the domain
    let mut result: Vec<f64> = vec![0 as f64; num_elements + 1];
    let spacing: f64 = (end - start).abs()/num_elements as f64;
    for (index, value) in result.iter_mut().enumerate() {
        *value = index as f64 * spacing;
    }
    return result;
}

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

    
}
