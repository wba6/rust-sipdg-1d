use util::linespace::linespace;

/// Generates a 1D discontinuous-Galerkin mesh on `[domain_a, domain_b]`.
///
/// This builds a uniform partition of the domain into `num_elements` elements, then
/// creates element-local degrees of freedom (DoFs) where each element owns its
/// left and right endpoint (DG nodes are not shared between neighboring elements).
///
/// Returns `(h_elem, x_dof)` where:
/// - `h_elem[i]` is the length of element `i`
/// - `x_dof` has length `2 * num_elements` and stores the DoF coordinates as
///   `[x_L0, x_R0, x_L1, x_R1, ..., x_L{N-1}, x_R{N-1}]`
pub fn generate_mesh(domain_a: f64, domain_b: f64, num_elements: usize) -> (Vec<f64>,Vec<f64>) {

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

    (h_elem, x_dof)
}
