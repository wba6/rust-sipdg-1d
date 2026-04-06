use util::linespace::linespace;
use crate::pde::BasisOrder;

/// Generates a 1D discontinuous-Galerkin mesh on `[domain_a, domain_b]`.
///
/// This builds a uniform partition of the domain into `num_elements` elements, then
/// creates element-local degrees of freedom (DoFs) based on the basis order.
///
/// Returns `(h_elem, x_dof)` where:
/// - `h_elem[i]` is the length of element `i`
/// - `x_dof` stores the DoF coordinates for each element sequentially.
pub fn generate_mesh(domain_a: f64, domain_b: f64, num_elements: usize, order: BasisOrder) -> (Vec<f64>, Vec<f64>) {
    // Generate evenly spaced points for the element interfaces
    let x_interface: Vec<f64> = linespace(domain_a, domain_b, num_elements + 1);

    let n_nodes = order.num_nodes();
    let n_dof: usize = n_nodes * num_elements;
    let mut x_dof: Vec<f64> = vec![0.0; n_dof];
    let mut h_elem: Vec<f64> = vec![0.0; num_elements];

    for i in 0..num_elements {
        let x_l: f64 = x_interface[i];
        let x_r: f64 = x_interface[i + 1];
        h_elem[i] = x_r - x_l;

        match order {
            BasisOrder::Linear => {
                x_dof[2 * i] = x_l;
                x_dof[2 * i + 1] = x_r;
            }
            BasisOrder::Quadratic => {
                x_dof[3 * i] = x_l;
                x_dof[3 * i + 1] = (x_l + x_r) / 2.0;
                x_dof[3 * i + 2] = x_r;
            }
        }
    }

    (h_elem, x_dof)
}
