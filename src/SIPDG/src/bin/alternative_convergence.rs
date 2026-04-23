use sipdg::{PdeProblem, SipdgAssembler, DirichletBC, NeumannBC, generate_mesh, pde::BasisOrder, pde::Side};
use util::cg::cg;
use util::quadrature::get_gauss_legendre_3pts;
use std::f64::consts::PI;

/// Problem: - (a p')' + q p = f on [0, 1]
/// a = 1, q = 1
/// Exact solution: p(x) = cos(pi * x)
/// f(x) = (pi^2 + 1) * cos(pi * x)
/// BCs:
/// x = 0: p(0) = 1 (Dirichlet)
/// x = 1: p'(1) = -pi * sin(pi) = 0 => a(1)p'(1) = 0 (Neumann)
struct CosineProblem;

impl PdeProblem for CosineProblem {
    fn a(&self, _x: f64) -> f64 { 1.0 }
    fn q(&self, _x: f64) -> f64 { 1.0 }
    fn f(&self, x: f64) -> f64 { (PI.powi(2) + 1.0) * (PI * x).cos() }
}

fn exact_soln(x: f64) -> f64 { (PI * x).cos() }

fn run_study(order: BasisOrder, sizes: &[usize], penalty: f64) -> Vec<(usize, f64)> {
    let mut results = Vec::new();

    for &num_elements in sizes {
        let (h_elem, x_dof) = generate_mesh(0.0, 1.0, num_elements, order);
        let mut assembler = SipdgAssembler::new(h_elem, x_dof, penalty, order);
        let prob = CosineProblem;

        assembler.assemble_volume(&prob);
        assembler.assemble_interfaces(&prob);

        let left_bc = DirichletBC { value: 1.0 };
        let right_bc = NeumannBC { value: 0.0 };
        
        let (mut mat, mut global_rhs) = assembler.assemble_to_global();
        assembler.apply_boundaries(&prob, &left_bc, &right_bc, &mut mat, &mut global_rhs);

        let op = assembler.matrix_free_op(&prob, &left_bc, &right_bc);
        let p = cg(&op, &global_rhs, 1e-12, 20000);

        let quad_pts = get_gauss_legendre_3pts();
        let n_nodes = order.num_nodes();

        let l2_err_sq: f64 = assembler.elements.iter().map(|elem| {
            let det_j = elem.jacobian();
            let mut local_err_sq = 0.0;
            for pt in &quad_pts {
                let x = elem.map_to_physical(pt.xi);
                let mut u_h = 0.0;
                for i in 0..n_nodes {
                    u_h += p[elem.nodes[i]] * elem.phi(i, pt.xi);
                }
                let u_ex = exact_soln(x);
                let err = u_h - u_ex;
                local_err_sq += pt.weight * err * err * det_j;
            }
            local_err_sq
        }).sum();

        let l2_error = l2_err_sq.sqrt();
        results.push((num_elements, l2_error));
        println!("N = {}, L2 Error = {:.6e}", num_elements, l2_error);
    }
    results
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sizes = [20, 40, 80, 160, 320];
    let penalty = 20.0;

    println!("--- Alternative Problem: Cosine with Dirichlet-Neumann BCs ---");
    println!("--- Linear Elements (p=1) ---");
    let results_p1 = run_study(BasisOrder::Linear, &sizes, penalty);
    
    println!("\nN, Error, Rate");
    for i in 0..results_p1.len() {
        if i == 0 {
            println!("{}, {:.6e}, -", results_p1[i].0, results_p1[i].1);
        } else {
            let rate = (results_p1[i-1].1 / results_p1[i].1).log2();
            println!("{}, {:.6e}, {:.2}", results_p1[i].0, results_p1[i].1, rate);
        }
    }

    println!("\n--- Quadratic Elements (p=2) ---");
    let results_p2 = run_study(BasisOrder::Quadratic, &sizes, penalty);

    println!("\nN, Error, Rate");
    for i in 0..results_p2.len() {
        if i == 0 {
            println!("{}, {:.6e}, -", results_p2[i].0, results_p2[i].1);
        } else {
            let rate = (results_p2[i-1].1 / results_p2[i].1).log2();
            println!("{}, {:.6e}, {:.2}", results_p2[i].0, results_p2[i].1, rate);
        }
    }

    Ok(())
}
