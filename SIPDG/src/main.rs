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


    // Penalty parameter for stability
    let penalty_param: u32 = 10;

    // ------------------- Generate Mesh --------------------

}
